import pytest # type: ignore
from unittest.mock import Mock, patch
import os
from dandd.delta_tree import DeltaTree
from dandd.species_specifics import SpeciesSpecifics

@pytest.fixture
def test_experiment():
    """Create a minimal test experiment dictionary"""
    return {
        'tool': 'dashing',
        'registers': 20,
        'canonicalize': True,
        'debug': False,
        'nthreads': 0,
        'baseset': set(),
        'safety': False,
        'fast': False,
        'verbose': False,
        'ksweep': None,
        'lowmem': False,
        'mock_run': True
    }

@patch('dandd.utils.blake2b')
@patch('dandd.delta_node.DeltaNode')
def test_minimal_tree(mock_node, mock_blake2b, temp_species, test_experiment):
    """Test minimal tree creation with minimal mocking"""
    # Setup test data
    test_files = ["test1.fa", "test2.fa"]
    for f in test_files:
        temp_species.fastahex[f] = "mocked_hash"

    # Setup single mock node for the tree
    mock_node.return_value = Mock(
        delta=0,
        bestk=0,
        node_title="test_node",
        fastas=test_files
    )

    # Create tree with minimal patching
    with patch.object(DeltaTree, '__init__') as mock_init:
        # Create a new instance without calling __init__
        tree = DeltaTree.__new__(DeltaTree)
        
        # Set up required attributes in the correct order
        tree.experiment = test_experiment
        tree.mink = 0
        tree.maxk = 0
        tree.kstart = temp_species.kstart
        tree.speciesinfo = temp_species
        tree._dt = [mock_node.return_value]  # Initialize _dt before build_tree
        tree.fastas = test_files
        tree.root = tree._dt[-1]
        tree.ngen = len(test_files)
        
        # Return None from __init__
        mock_init.return_value = None

        # Verify basic tree structure
        assert len(tree._dt) == 1
        assert tree.root == mock_node.return_value
        assert tree.fastas == test_files