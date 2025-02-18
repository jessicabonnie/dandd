import pytest # type: ignore
import os
from unittest.mock import Mock, patch
from dandd.delta_tree_node import DeltaTreeNode
from dandd.sketch_filepath import SketchFilePath

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
def test_find_delta(mock_run, temp_species, test_experiment):
    """Test find_delta method"""
    # Set up mock subprocess behavior
    mock_run.return_value.stdout = "#Path\tSize (est.)\n/path/to/sketch\t1000.0\n"
    
    node = DeltaTreeNode("test_node", [], temp_species, test_experiment)
    node.find_delta(10)
    
    assert node.bestk > 0
    assert node.delta > 0
    assert hasattr(node, 'card')

@patch('subprocess.run')
def test_node_ksweep(mock_run, temp_species, test_experiment):
    """Test node_ksweep method"""
    # Set up mock subprocess behavior
    mock_run.return_value.stdout = "#Path\tSize (est.)\n/path/to/sketch\t1000.0\n"
    
    node = DeltaTreeNode("test_node", [], temp_species, test_experiment)
    node.node_ksweep(5, 7)
    
    # Check that sketches were created for each k
    for k in range(5, 8):
        assert node.ksketches[k] is not None
        assert node.ksketches[k].kval == k
    
    assert node.mink == 5
    assert node.maxk == 7 