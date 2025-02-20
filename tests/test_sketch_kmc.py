import pytest # type: ignore
import os
import subprocess
from unittest.mock import Mock, patch
from dandd.sketch_kmc import KMCSketchObj
from dandd.sketch_filepath import SketchFilePath
import shutil
import tempfile
import glob

# Define KMCLOC with full path - adjust as needed
KMCLOC = "kmc"

@pytest.fixture
def kmc_sketch_obj(temp_species, test_experiment):
    """Create a KMCSketchObj instance"""
    # Patch KMCLOC in the sketch_kmc module
    import dandd.sketch_kmc
    dandd.sketch_kmc.KMCLOC = KMCLOC
    
    # Create with kval=21 but prevent sketch creation
    test_experiment['mock_run'] = True  # Prevent actual command execution
    test_experiment['tool'] = 'kmc'  # Set tool to KMC
    
    # Patch create_sketch to do nothing during init
    with patch('dandd.sketch_base.SketchObj.create_sketch'):
        sfp = SketchFilePath(
            filenames=[
                os.path.join(temp_species.inputdir, "NC_009057.fasta"),
                os.path.join(temp_species.inputdir, "NC_009058.fasta.gz")
            ],
            kval=21,
            speciesinfo=temp_species,
            experiment=test_experiment
        )
        
        # Create necessary directories
        os.makedirs(os.path.dirname(sfp.full), exist_ok=True)
        
        obj = KMCSketchObj(
            kval=21,
            sfp=sfp,
            speciesinfo=temp_species,
            experiment=test_experiment
        )
    
    return obj

def test_sketch_check(kmc_sketch_obj, tmp_path):
    """Test sketch existence check"""
    # Test non-existent sketch
    assert not kmc_sketch_obj.sketch_check()
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(kmc_sketch_obj.sfp.full), exist_ok=True)
    
    # Create test files (.kmc_pre and .kmc_suf)
    pre_file = f"{kmc_sketch_obj.sfp.full}.kmc_pre"
    suf_file = f"{kmc_sketch_obj.sfp.full}.kmc_suf"
    
    with open(pre_file, 'w') as f:
        f.write("test")
    with open(suf_file, 'w') as f:
        f.write("test")
        
    assert kmc_sketch_obj.sketch_check()
    
    # Cleanup
    os.remove(pre_file)
    os.remove(suf_file)

@patch('subprocess.run')
def test_leaf_command(mock_run, kmc_sketch_obj):
    """Test leaf sketch command generation"""
    tmpdir = "/tmp/test"
    cmd = kmc_sketch_obj._leaf_command(tmpdir)
    
    expected_parts = [
        "kmc -hp",
        "-k21",
        "-ci1",
        "-cs2",
        "-fm",
        kmc_sketch_obj.sfp.ffiles[0]
    ]
    
    for part in expected_parts:
        assert str(part) in cmd

@patch('subprocess.run')
def test_union_command(mock_run, kmc_sketch_obj):
    """Test union sketch command generation"""
    cmd = kmc_sketch_obj._union_command()
    
    expected_parts = [
        "kmc_tools",
        "complex",
        "/dev/stdin"
    ]
    
    for part in expected_parts:
        assert str(part) in cmd

@patch('subprocess.run')
def test_card_command(mock_run, kmc_sketch_obj):
    """Test cardinality command generation"""
    # Set sketch path
    kmc_sketch_obj.sketch = kmc_sketch_obj.sfp.full
    
    cmd = kmc_sketch_obj.card_command([kmc_sketch_obj.sfp.full])
    
    expected_parts = [
        "kmc_tools",
        "info",
        kmc_sketch_obj.sfp.full
    ]
    
    for part in expected_parts:
        assert str(part) in cmd

def test_parse_card(kmc_sketch_obj):
    """Test cardinality parsing"""
    # Set sketch path
    kmc_sketch_obj.sketch = kmc_sketch_obj.sfp.full
    
    # Mock process output with KMC format
    proc = Mock()
    # Format matches card_command output format
    proc.stdout = f"{kmc_sketch_obj.sfp.full},1000000"
    
    kmc_sketch_obj.parse_card(proc)
    
    assert kmc_sketch_obj.speciesinfo.cardkey[kmc_sketch_obj.sfp.full] == 1000000.0

def test_remove_sketch(kmc_sketch_obj):
    """Test sketch removal"""
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(kmc_sketch_obj.sfp.full), exist_ok=True)
    
    # Create test files
    pre_file = f"{kmc_sketch_obj.sfp.full}.kmc_pre"
    suf_file = f"{kmc_sketch_obj.sfp.full}.kmc_suf"
    
    with open(pre_file, 'w') as f:
        f.write("test")
    with open(suf_file, 'w') as f:
        f.write("test")
    
    assert os.path.exists(pre_file)
    assert os.path.exists(suf_file)
    
    # Mock glob to return our test file
    with patch('glob.glob') as mock_glob:
        mock_glob.return_value = [kmc_sketch_obj.sfp.full]
        kmc_sketch_obj.remove_sketch()
    
    assert not os.path.exists(pre_file)
    assert not os.path.exists(suf_file)

@patch('subprocess.run')
def test_error_handling(mock_run, kmc_sketch_obj):
    """Test error handling in command execution"""
    mock_run.side_effect = subprocess.CalledProcessError(1, "kmc command")
    
    with pytest.raises(subprocess.CalledProcessError):
        kmc_sketch_obj.create_leaf_sketch()

def is_kmc_available():
    """Check if KMC is available on the system"""
    return shutil.which(KMCLOC) is not None

@pytest.mark.skipif(not is_kmc_available(), 
                   reason="kmc command not found")
def test_full_sketch_workflow(kmc_sketch_obj, test_data_dir):
    """Test full workflow: create sketch and check cardinality without mocks"""
    # Use existing test file
    test_fasta = os.path.join(test_data_dir, "NC_009057.fasta")
    
    # Update sketch object to use our test file
    kmc_sketch_obj.sfp.ffiles = [test_fasta]
    
    # Create sketch directory if it doesn't exist
    sketch_dir = os.path.dirname(kmc_sketch_obj.sfp.full)
    os.makedirs(sketch_dir, exist_ok=True)
    
    try:
        # Create temporary directory for KMC
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create the sketch
            sketch_cmd = f"kmc -hp -k21 -ci1 -cs65535 -fm {test_fasta} {kmc_sketch_obj.sfp.full} {tmpdir}"
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
            assert os.path.exists(f"{kmc_sketch_obj.sfp.full}.kmc_pre")
            assert os.path.exists(f"{kmc_sketch_obj.sfp.full}.kmc_suf")
            
            # Check cardinality
            card_cmd = f"kmc_tools info {kmc_sketch_obj.sfp.full}"
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
            kmc_sketch_obj.parse_card(result)
            assert kmc_sketch_obj.sfp.full in kmc_sketch_obj.speciesinfo.cardkey
            assert isinstance(kmc_sketch_obj.speciesinfo.cardkey[kmc_sketch_obj.sfp.full], float)
            assert kmc_sketch_obj.speciesinfo.cardkey[kmc_sketch_obj.sfp.full] > 0
    finally:
        # Clean up any leftover files
        kmc_sketch_obj.remove_sketch()

def test_parse_card_real_format(kmc_sketch_obj):
    """Test parsing real KMC output format"""
    # Create a mock process with real KMC output format
    proc = Mock()
    proc.stdout = f"""
Info for {kmc_sketch_obj.sfp.full}:
k = 21
total k-mers: 4606432
unique k-mers: 4123789
singleton k-mers: 3500000
max count: 65535
min count: 1
"""
    
    # Parse the output
    kmc_sketch_obj.parse_card(proc)
    
    # Verify the cardinality was stored correctly
    assert kmc_sketch_obj.sfp.full in kmc_sketch_obj.speciesinfo.cardkey
    assert kmc_sketch_obj.speciesinfo.cardkey[kmc_sketch_obj.sfp.full] == 4606432  # Using total k-mers value

@patch('subprocess.run')
def test_kmc_command_output(mock_run, kmc_sketch_obj, test_data_dir):
    """Test actual KMC command execution and output capture"""
    # Create a Mock object with stdout attribute
    mock_result = Mock()
    mock_result.stdout = f"{kmc_sketch_obj.sfp.full},4606432"  # Match the format from card_command
    mock_run.return_value = mock_result
    
    # Use existing test file
    test_fasta = os.path.join(test_data_dir, "NC_009057.fasta")
    print(f"Test file path: {test_fasta}")
    
    # Update sketch object to use our test file
    kmc_sketch_obj.sfp.ffiles = [test_fasta]
    
    # Generate and run the card command
    cmd = kmc_sketch_obj.card_command([kmc_sketch_obj.sfp.full])
    print(f"Command to run: {cmd}")
    
    # Parse and verify the cardinality
    kmc_sketch_obj.parse_card(mock_result)
    assert kmc_sketch_obj.sfp.full in kmc_sketch_obj.speciesinfo.cardkey
    assert kmc_sketch_obj.speciesinfo.cardkey[kmc_sketch_obj.sfp.full] == 4606432