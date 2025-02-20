import pytest # type: ignore
import os
from unittest.mock import Mock, patch
from dandd.delta_node import DeltaTreeNode
# from dandd.sketch_filepath import SketchFilePath

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
    """Test initialization of DeltaTreeNode"""
    node = DeltaTreeNode(
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
    node = DeltaTreeNode("test_node", [], temp_species, test_experiment)
    node.assign_progeny()
    
    assert node.progeny == [node]
    assert node.fastas == ["test_node"]
    assert node.node_title == "test_node"  # Since no extension in test name

def test_lt_comparison(temp_species, test_experiment):
    """Test less than comparison between nodes"""
    node1 = DeltaTreeNode("node1", [], temp_species, test_experiment)
    node2 = DeltaTreeNode("node2", [], temp_species, test_experiment)
    node2.ngen = 2
    
    assert node1 < node2
    assert not node2 < node1

def test_find_delta_helper_dashing_limit(temp_species, test_experiment):
    """Test that find_delta_helper raises error for dashing with k > 32"""
    node = DeltaTreeNode("test_node", [], temp_species, test_experiment)
    
    with pytest.raises(ValueError) as exc_info:
        node.find_delta_helper(33)
    assert "Exploratory k value is too high for dashing" in str(exc_info.value)

def test_node_with_children(temp_species, test_experiment):
    """Test node with children initialization"""
    child1 = DeltaTreeNode("child1", [], temp_species, test_experiment)
    child2 = DeltaTreeNode("child2", [], temp_species, test_experiment)
    parent = DeltaTreeNode("parent", [child1, child2], temp_species, test_experiment)
    
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
    
    node = DeltaTreeNode("test_node", [], temp_species, experiment)
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
    node = DeltaTreeNode(test_file, [test_file], temp_species, test_experiment)
    
    # Start at k=20 to avoid hitting dashing's k=32 limit during exploration
    # Mock the recursive calls to find_delta_helper
    with patch.object(DeltaTreeNode, 'find_delta_helper') as mock_helper:
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
    
    # Create node with real file path
    node = DeltaTreeNode(test_file, [test_file], temp_species, test_experiment)
    
    # Use k values within dashing's limit
    node.node_ksweep(5, 30)
    
    # Add assertions
    assert len(node.ksketches) > 0
    assert all(sketch.kval <= 32 for sketch in node.ksketches if sketch is not None)  # Verify all k values are within limit 