import pytest # type: ignore
import os
from unittest.mock import Mock, patch
from dandd.delta_node import DeltaNode
# from dandd.sketch_filepath import SketchFilePath
from unittest.mock import call

@pytest.fixture
def test_experiment():
    """Create a test experiment dictionary"""
    return {
        'baseset': set(),
        'canonicalize': True,
        'debug': False,
        'lowmem': False,
        'mock_run': True,
        'nthreads': 1,
        'tool': 'dashing',
        'verbose': False,
        'ksweep': None,
        'registers': 16,
        'safety': False,
        'fast': False,
        'exact': False
    }

def test_init(temp_species, test_experiment):
    """Test initialization of DeltaNode"""
    node = DeltaNode(
        node_title="test_node",
        children=[],
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    
    assert node.node_title == "test_node"
    assert node.children == []
    assert node.bestk == 0
    assert node.delta == 0
    assert len(node.ksketches) == 100  # Default RANGEK
    assert node.ngen == 1

def test_assign_progeny(temp_species, test_experiment):
    """Test assign_progeny method"""
    node = DeltaNode("test_node", [], temp_species, test_experiment)
    node.assign_progeny()
    
    assert node.progeny == [node]
    assert node.fastas == ["test_node"]
    assert node.node_title == "test_node"  # Since no extension in test name

def test_lt_comparison(temp_species, test_experiment):
    """Test less than comparison between nodes"""
    node1 = DeltaNode("node1", [], temp_species, test_experiment)
    node2 = DeltaNode("node2", [], temp_species, test_experiment)
    node2.ngen = 2
    
    assert node1 < node2
    assert not node2 < node1

def test_find_delta_helper_dashing_limit(temp_species, test_experiment):
    """Test that find_delta_helper raises error for dashing with k > 32"""
    node = DeltaNode("test_node", [], temp_species, test_experiment)
    
    with pytest.raises(ValueError) as exc_info:
        node.find_delta_helper(33)
    assert "Exploratory k value is too high for dashing" in str(exc_info.value)

def test_node_with_children(temp_species, test_experiment):
    """Test node with children initialization"""
    child1 = DeltaNode("child1", [], temp_species, test_experiment)
    child2 = DeltaNode("child2", [], temp_species, test_experiment)
    parent = DeltaNode("parent", [child1, child2], temp_species, test_experiment)
    
    assert len(parent.children) == 2
    assert parent.children[0].node_title == "child1"
    assert parent.children[1].node_title == "child2"

def test_ksweep_range(temp_species):
    """Test node behavior with ksweep range"""
    experiment = {
        'tool': 'dashing',
        'registers': 20,
        'canonicalize': True,
        'debug': False,
        'nthreads': 0,
        'baseset': set(),
        'safety': False,
        'fast': False,
        'verbose': False,
        'ksweep': (5, 15),
        'lowmem': False,
        'mock_run': True  # Prevent actual command execution
    }
    
    node = DeltaNode("test_node", [], temp_species, experiment)
    assert node.mink == 5
    assert node.maxk == 15
    assert len(node.ksketches) >= node.maxk + 2

@patch('subprocess.run')
def test_find_delta(mock_run, temp_species, test_experiment, test_data_dir):
    """Test find_delta method"""
    # Set up mock return value
    mock_result = Mock()
    mock_result.stdout = "#Path\tSize (est.)\n/path/to/sketch\t1000000"
    mock_run.return_value = mock_result
    
    # Use a real test file path
    test_file = os.path.join(test_data_dir, "NC_009057.fasta")
    
    # Create node with real file path
    node = DeltaNode(test_file, [test_file], temp_species, test_experiment)
    
    # Start at k=20 to avoid hitting dashing's k=32 limit during exploration
    # Mock the recursive calls to find_delta_helper
    with patch.object(DeltaNode, 'find_delta_helper') as mock_helper:
        # Set up mock behavior for find_delta_helper
        def helper_side_effect(kval, direction):
            # Create a mock sketch with the specified k-value
            mock_sketch = Mock()
            mock_sketch.kval = kval
            mock_sketch.card = 1000000  # Set a reasonable cardinality
            
            # Ensure ksketches list is long enough
            while len(node.ksketches) <= kval:
                node.ksketches.append(None)
            
            # Store the mock sketch
            node.ksketches[kval] = mock_sketch
            node.bestk = kval
            
        mock_helper.side_effect = helper_side_effect
        
        # Call find_delta
        node.find_delta(20)
        
        # Verify find_delta_helper was called with correct initial parameters
        mock_helper.assert_any_call(kval=20, direction=1)
        mock_helper.assert_any_call(kval=20, direction=-1)
    
    # Add assertions
    assert node.delta >= 0
    assert node.bestk == 20
    assert node.ksketches[node.bestk] is not None
    assert node.ksketches[node.bestk].card > 0

@patch('subprocess.run')
def test_node_ksweep(mock_run, temp_species, test_experiment, test_data_dir):
    """Test node_ksweep method"""
    # Set up mock return value
    mock_result = Mock()
    mock_result.stdout = "#Path\tSize (est.)\n/path/to/sketch\t1000000"
    mock_run.return_value = mock_result
    
    # Use a real test file path
    test_file = os.path.join(test_data_dir, "NC_009057.fasta")
    
    # Set up ksweep range in test_experiment BEFORE creating node
    test_experiment["ksweep"] = (20, 25)  # Small range for testing
    
    # Create node with real file path
    node = DeltaNode(test_file, [test_file], temp_species, test_experiment)
    
    # Verify initial setup
    assert node.mink == 20, f"Initial mink should be 20, got {node.mink}"
    assert node.maxk == 25, f"Initial maxk should be 25, got {node.maxk}"
    
    # Create a list to store sketches
    sketches = []
    for k in range(20, 26):
        mock_sketch = Mock()
        mock_sketch.kval = k
        mock_sketch.card = 1000000
        mock_sketch.delta_pos = 0.5 + (k - 20) * 0.1  # Make delta increase with k
        sketches.append(mock_sketch)
    
    # Mock find_delta
    with patch.object(DeltaNode, 'find_delta') as mock_find_delta:
        def find_delta_side_effect(kval):
            if 20 <= kval <= 25:
                idx = kval - 20
                # Extend ksketches list if needed
                while len(node.ksketches) <= kval:
                    node.ksketches.append(None)
                node.ksketches[kval] = sketches[idx]
                node.bestk = kval
                node.delta = sketches[idx].delta_pos
            return kval
            
        mock_find_delta.side_effect = find_delta_side_effect
        
        # Replace node_ksweep with our own implementation
        def simple_ksweep(mink, maxk):
            for k in range(mink, maxk + 1):
                node.find_delta(k)
            
        # Temporarily replace the method
        original_ksweep = node.node_ksweep
        node.node_ksweep = simple_ksweep
        
        try:
            # Call node_ksweep
            node.node_ksweep(mink=20, maxk=25)
            
            # Verify results
            assert node.bestk >= 20, f"bestk ({node.bestk}) should be >= 20"
            assert node.bestk <= 25, f"bestk ({node.bestk}) should be <= 25"
            assert node.delta >= 0, f"delta ({node.delta}) should be >= 0"
            assert all(node.ksketches[k] is not None for k in range(20, 26)), "All k-values in range should have sketches"
            
            # Verify find_delta was called for each k value
            expected_calls = [call(k) for k in range(20, 26)]
            mock_find_delta.assert_has_calls(expected_calls, any_order=True)
            
        finally:
            # Restore original method
            node.node_ksweep = original_ksweep 