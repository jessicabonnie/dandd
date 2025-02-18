import pytest # type: ignore
import os
import subprocess
from unittest.mock import Mock, patch
from dandd.sketch_dashing import DashSketchObj
from dandd.sketch_filepath import SketchFilePath
import glob
import shutil

# Define DASHINGLOC with full path
DASHINGLOC = "/home/jbonnie1/lib/dashing/dashing"

@pytest.fixture
def dash_sketch_obj(temp_species, test_experiment):
    """Create a DashSketchObj instance"""
    # Patch DASHINGLOC in the sketch_dashing module
    import dandd.sketch_dashing
    dandd.sketch_dashing.DASHINGLOC = DASHINGLOC
    
    # Create with kval=21 but prevent sketch creation
    test_experiment['mock_run'] = True  # Prevent actual command execution
    
    # Patch create_sketch to do nothing during init
    with patch('dandd.sketch_base.SketchObj.create_sketch'):
        sfp = SketchFilePath(
            filenames=[
                os.path.join(temp_species.inputdir, "NC_009057.fasta"),
                os.path.join(temp_species.inputdir, "NC_009058.fasta.gz")
            ],
            kval=21,  # Use actual kval
            speciesinfo=temp_species,
            experiment=test_experiment
        )
        
        # Create necessary directories
        os.makedirs(os.path.dirname(sfp.full), exist_ok=True)
        
        obj = DashSketchObj(
            kval=21,  # Use actual kval
            sfp=sfp,
            speciesinfo=temp_species,
            experiment=test_experiment
        )
    
    return obj

def test_sketch_check(dash_sketch_obj, tmp_path):
    """Test sketch existence check"""
    # Test non-existent sketch
    assert not dash_sketch_obj.sketch_check()
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(dash_sketch_obj.sfp.full), exist_ok=True)
    
    # Test existing sketch
    with open(dash_sketch_obj.sfp.full, 'w') as f:
        f.write("test")
    assert dash_sketch_obj.sketch_check()
    os.remove(dash_sketch_obj.sfp.full)

@patch('subprocess.run')
def test_leaf_command(mock_run, dash_sketch_obj):
    """Test leaf sketch command generation"""
    tmpdir = "/tmp/test"
    cmd = dash_sketch_obj._leaf_command(tmpdir)
    
    expected_parts = [
        DASHINGLOC,
        "sketch",
        "-k21",  # Changed format
        "-S", "16",  # Changed from -p to -S
        dash_sketch_obj.sfp.ffiles[0],  # input file
    ]
    
    for part in expected_parts:
        assert str(part) in cmd

@patch('subprocess.run')
def test_union_command(mock_run, dash_sketch_obj):
    """Test union sketch command generation"""
    cmd = dash_sketch_obj._union_command()
    
    expected_parts = [
        DASHINGLOC,
        "union",
        "-z",  # Added actual flag
        "-o", dash_sketch_obj.sfp.full  # output file
    ]
    
    for part in expected_parts:
        assert str(part) in cmd

@patch('subprocess.run')
def test_card_command(mock_run, dash_sketch_obj):
    """Test cardinality command generation"""
    # Set sketch path
    dash_sketch_obj.sketch = dash_sketch_obj.sfp.full
    
    cmd = dash_sketch_obj.card_command()
    
    expected_parts = [
        DASHINGLOC,
        "card",
        "--presketched",
        dash_sketch_obj.sfp.full
    ]
    
    for part in expected_parts:
        assert str(part) in cmd

def test_parse_card(dash_sketch_obj):
    """Test cardinality parsing"""
    # Set sketch path
    dash_sketch_obj.sketch = dash_sketch_obj.sfp.full
    
    # Mock process output with proper TSV format
    proc = Mock()
    proc.stdout = "#Path\tSize (est.)\n" + dash_sketch_obj.sfp.full + "\t1234.56\n"
    
    dash_sketch_obj.parse_card(proc)
    
    assert dash_sketch_obj.speciesinfo.cardkey[dash_sketch_obj.sfp.full] == 1234.56

@patch('glob.glob', autospec=True)  # Use autospec and patch glob.glob directly
def test_remove_sketch(mock_glob, dash_sketch_obj):
    """Test sketch removal"""
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(dash_sketch_obj.sfp.full), exist_ok=True)
    
    # Create test file
    with open(dash_sketch_obj.sfp.full, 'w') as f:
        f.write("test")
    
    assert os.path.exists(dash_sketch_obj.sfp.full)
    
    # Set up mock return value
    mock_glob.return_value = [dash_sketch_obj.sfp.full]
    
    # Import and patch glob in sketch_dashing
    import dandd.sketch_dashing
    dandd.sketch_dashing.glob = glob
    
    dash_sketch_obj.remove_sketch()
    
    assert not os.path.exists(dash_sketch_obj.sfp.full)

@patch('subprocess.run')
def test_error_handling(mock_run, dash_sketch_obj):
    """Test error handling in command execution"""
    mock_run.side_effect = subprocess.CalledProcessError(1, "dashing command")
    
    with pytest.raises(subprocess.CalledProcessError):
        dash_sketch_obj.create_leaf_sketch()

def is_dashing_available():
    """Check if dashing is available on the system"""
    return shutil.which(DASHINGLOC) is not None

@pytest.fixture
def test_data_dir():
    """Path to test data directory"""
    return os.path.join(os.path.dirname(__file__), 'data')

@pytest.mark.skipif(not is_dashing_available(), 
                   reason="dashing command not found")
def test_full_sketch_workflow(dash_sketch_obj, test_data_dir):
    """Test full workflow: create sketch and check cardinality without mocks"""
    # Use existing test file
    test_fasta = os.path.join(test_data_dir, "NC_009057.fasta")
    
    # Update sketch object to use our test file
    dash_sketch_obj.sfp.ffiles = [test_fasta]
    
    # Create sketch directory if it doesn't exist
    sketch_dir = os.path.dirname(dash_sketch_obj.sfp.full)
    os.makedirs(sketch_dir, exist_ok=True)
    
    try:
        # Create the sketch
        sketch_cmd = f"{DASHINGLOC} sketch -k21 -S16 {test_fasta} -o {dash_sketch_obj.sfp.full}"
        print(f"Sketch command: {sketch_cmd}")
        
        result = subprocess.run(
            sketch_cmd,
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )
        print(f"Sketch output: {result.stdout}")
        
        # Verify sketch was created
        assert os.path.exists(dash_sketch_obj.sfp.full)
        
        # Check cardinality
        card_cmd = f"{DASHINGLOC} card --presketched {dash_sketch_obj.sfp.full}"
        print(f"Card command: {card_cmd}")
        
        result = subprocess.run(
            card_cmd,
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )
        print(f"Card output: {result.stdout}")
        
        # Parse and verify cardinality
        dash_sketch_obj.parse_card(result)
        assert dash_sketch_obj.sfp.full in dash_sketch_obj.speciesinfo.cardkey
        assert isinstance(dash_sketch_obj.speciesinfo.cardkey[dash_sketch_obj.sfp.full], float)
        assert dash_sketch_obj.speciesinfo.cardkey[dash_sketch_obj.sfp.full] > 0
        
    finally:
        # Clean up
        if os.path.exists(dash_sketch_obj.sfp.full):
            os.remove(dash_sketch_obj.sfp.full)

@pytest.mark.skipif(not is_dashing_available(), 
                   reason="dashing command not found")
def test_real_dashing_command(dash_sketch_obj, test_data_dir):
    """Test actual dashing command execution without mocks"""
    # Use existing test file
    test_fasta = os.path.join(test_data_dir, "NC_009057.fasta")
    
    # Run dashing directly to get cardinality
    cmd = f"{DASHINGLOC} card {test_fasta}"
    print(f"Command to run: {cmd}")
    
    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        capture_output=True,
        text=True
    )
    print(f"Command output: {result.stdout}")
    
    # Verify output format matches actual dashing output
    assert "#Path" in result.stdout
    assert "Size (est.)" in result.stdout
    
    # Parse the output
    dash_sketch_obj.parse_card(result)
    assert test_fasta in dash_sketch_obj.speciesinfo.cardkey
    assert isinstance(dash_sketch_obj.speciesinfo.cardkey[test_fasta], float)
    assert dash_sketch_obj.speciesinfo.cardkey[test_fasta] > 0

@pytest.mark.skipif(not is_dashing_available(), 
                   reason="dashing command not found")
def test_sketch_and_card(dash_sketch_obj, test_data_dir):
    """Test creating a sketch and checking its cardinality"""
    # Use existing test file
    test_fasta = os.path.join(test_data_dir, "NC_009057.fasta")
    
    # Update sketch object to use our test file
    dash_sketch_obj.sfp.ffiles = [test_fasta]
    
    # Create sketch directory if it doesn't exist
    os.makedirs(os.path.dirname(dash_sketch_obj.sfp.full), exist_ok=True)
    
    # Generate and run the sketch command
    sketch_cmd = f"{DASHINGLOC} sketch -k21 -S16 {test_fasta} -o {dash_sketch_obj.sfp.full}"
    print(f"Sketch command: {sketch_cmd}")
    
    # Create the sketch
    with patch('subprocess.run', wraps=subprocess.run) as real_run:
        with patch.dict(dash_sketch_obj.experiment, {'mock_run': False}):
            result = real_run(
                sketch_cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            print(f"Sketch output: {result.stdout}")
    
    # Verify sketch was created
    assert os.path.exists(dash_sketch_obj.sfp.full)
    
    # Now check cardinality
    card_cmd = dash_sketch_obj.card_command([dash_sketch_obj.sfp.full])
    print(f"Card command: {card_cmd}")
    
    with patch('subprocess.run', wraps=subprocess.run) as real_run:
        with patch.dict(dash_sketch_obj.experiment, {'mock_run': False}):
            result = real_run(
                card_cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            print(f"Card output: {result.stdout}")
    
    # Parse and verify cardinality
    dash_sketch_obj.parse_card(result)
    assert dash_sketch_obj.sfp.full in dash_sketch_obj.speciesinfo.cardkey
    assert isinstance(dash_sketch_obj.speciesinfo.cardkey[dash_sketch_obj.sfp.full], float)
    assert dash_sketch_obj.speciesinfo.cardkey[dash_sketch_obj.sfp.full] > 0
    
    # Clean up
    os.remove(dash_sketch_obj.sfp.full)

def test_parse_card_real_format(dash_sketch_obj, test_data_dir):
    """Test parsing real dashing output format"""
    # Create a mock process with real dashing output format
    proc = Mock()
    proc.stdout = """#Path	Size (est.)	Card. Type	Exact?
/home/jbonnie1/dandd/dandd/tests/data/NC_009057.fasta	4606432	Hyperlogs	No
/home/jbonnie1/dandd/dandd/tests/data/NC_009058.fasta.gz	4123789	Hyperlogs	No"""

    # Parse the output
    dash_sketch_obj.parse_card(proc)

    # Verify the parsed values
    test_files = [
        os.path.join(test_data_dir, "NC_009057.fasta"),
        os.path.join(test_data_dir, "NC_009058.fasta.gz")
    ]
    
    # Check that both files were parsed
    for test_file in test_files:
        assert test_file in dash_sketch_obj.speciesinfo.cardkey
        assert isinstance(dash_sketch_obj.speciesinfo.cardkey[test_file], float)
    
    # Check specific values
    assert dash_sketch_obj.speciesinfo.cardkey[test_files[0]] == 4606432.0
    assert dash_sketch_obj.speciesinfo.cardkey[test_files[1]] == 4123789.0 

@patch('subprocess.run')
def test_dashing_command_output(mock_run, dash_sketch_obj, test_data_dir):
    """Test actual dashing command execution and output capture"""
    # Create a Mock object with stdout attribute
    mock_result = Mock()
    mock_result.stdout = """#Path\tSize (est.)\tCard. Type\tExact?
/home/jbonnie1/dandd/dandd/tests/data/NC_009057.fasta\t4606432\tHyperlogs\tNo"""
    mock_run.return_value = mock_result
    
    # Use existing test file
    test_fasta = os.path.join(test_data_dir, "NC_009057.fasta")
    print(f"Test file path: {test_fasta}")
    
    # Update sketch object to use our test file
    dash_sketch_obj.sfp.ffiles = [test_fasta]
    
    # Generate and run the card command
    cmd = dash_sketch_obj.card_command([test_fasta])
    print(f"Command to run: {cmd}")
    
    # Run command and capture output
    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        capture_output=True,
        text=True
    )
    
    # Verify output format
    assert "#Path" in mock_result.stdout
    assert "Size (est.)" in mock_result.stdout
    
    # Parse and verify the cardinality
    dash_sketch_obj.parse_card(mock_result)
    assert test_fasta in dash_sketch_obj.speciesinfo.cardkey
    assert isinstance(dash_sketch_obj.speciesinfo.cardkey[test_fasta], float)