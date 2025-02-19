import pytest # type: ignore
import subprocess
import tempfile
import os
from unittest.mock import Mock, patch
from dandd.sketch_base import SketchObj
from dandd.sketch_filepath import SketchFilePath
import hashlib  # Add this import

class MockSketchObj(SketchObj):
    """Mock implementation of SketchObj for testing"""
    _mock_sketch_check = None
    _mock_leaf_command = None
    _mock_remove_sketch = None
    _mock_card_command = None
    _mock_parse_card = None
    _mock_union_command = None
    
    def sketch_check(self, path=None) -> bool:
        return self._mock_sketch_check(path)

    def _leaf_command(self, tmpdir) -> str:
        return self._mock_leaf_command(tmpdir)

    def remove_sketch(self):
        return self._mock_remove_sketch()

    def card_command(self, sketchlist) -> str:
        return self._mock_card_command(sketchlist)

    def parse_card(self, proc):
        return self._mock_parse_card(proc)

    def _union_command(self) -> str:
        return self._mock_union_command()

@pytest.fixture
def mock_run():
    """Mock subprocess.run"""
    with patch('subprocess.run') as mock:
        mock.return_value = Mock(returncode=0)
        yield mock

@pytest.fixture
def mock_sketch_obj(temp_species, test_experiment, mock_run):
    """Create a MockSketchObj instance with mocked methods"""
    # Set up mock methods as class attributes first
    MockSketchObj._mock_sketch_check = Mock(return_value=False)
    MockSketchObj._mock_leaf_command = Mock(return_value="mock leaf command")
    MockSketchObj._mock_remove_sketch = Mock()
    MockSketchObj._mock_card_command = Mock(return_value="mock card command")
    MockSketchObj._mock_parse_card = Mock()
    MockSketchObj._mock_union_command = Mock(return_value="mock union command")
    
    # Create test files and add their hashes to fastahex
    test_files = [
        os.path.join(temp_species.inputdir, "test1.fa"),
        os.path.join(temp_species.inputdir, "test2.fasta.gz")
    ]
    
    # Create mock hashes for test files without trying to read them
    temp_species.fastahex = {}
    for test_file in test_files:
        # Create a deterministic mock hash using the filename
        mock_hash = hashlib.blake2b(os.path.basename(test_file).encode()).hexdigest()
        temp_species.fastahex[os.path.basename(test_file)] = mock_hash
    
    sfp = SketchFilePath(
        filenames=test_files,
        kval=0,  # Set kval to 0 to skip initialization
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    
    # Don't initialize cardkey here - let tests handle it
    temp_species.cardkey = {}
    
    obj = MockSketchObj(
        kval=0,  # Set kval to 0 to skip initialization
        sfp=sfp,
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    
    # Now set kval to 21 for testing
    obj.kval = 21
    
    return obj

def test_init(mock_sketch_obj):
    """Test initialization of SketchObj"""
    assert mock_sketch_obj.kval == 21
    assert mock_sketch_obj.sketch is None
    assert mock_sketch_obj.cmd is None
    assert mock_sketch_obj.delta_pos == 0
    assert mock_sketch_obj.card == 0
    assert isinstance(mock_sketch_obj.sfp, SketchFilePath)

@patch('subprocess.run')
def test_create_leaf_sketch(mock_run, mock_sketch_obj):
    """Test leaf sketch creation"""
    mock_run.return_value = Mock(returncode=0)
    
    mock_sketch_obj.create_leaf_sketch()
    
    assert mock_sketch_obj.cmd == "mock leaf command"
    mock_run.assert_called_once()
    assert mock_run.call_args[1]['shell'] is True
    assert mock_run.call_args[1]['check'] is True

@patch('subprocess.run')
def test_create_union_sketch(mock_run, mock_sketch_obj):
    """Test union sketch creation"""
    mock_run.return_value = Mock(returncode=0)
    
    mock_sketch_obj.create_union_sketch()
    
    assert mock_sketch_obj.cmd == "mock union command"
    mock_run.assert_called_once()
    assert mock_run.call_args[1]['shell'] is True
    assert mock_run.call_args[1]['check'] is True

@patch('subprocess.run')
def test_individual_card(mock_run, mock_sketch_obj):
    """Test individual cardinality calculation"""
    mock_run.return_value = Mock(
        returncode=0,
        stdout="mock output",
        stderr="",
    )
    
    mock_sketch_obj.individual_card()
    
    assert mock_sketch_obj._mock_card_command.called
    assert mock_sketch_obj._mock_parse_card.called
    mock_run.assert_called_once()

def test_ngen_property(mock_sketch_obj):
    """Test ngen property returns correct number of input files"""
    assert mock_sketch_obj.ngen == mock_sketch_obj.sfp.ngen

@patch('subprocess.run')
def test_subprocess_error_handling(mock_run, mock_sketch_obj):
    """Test error handling for subprocess failures"""
    mock_run.side_effect = subprocess.CalledProcessError(1, "mock command")
    
    with pytest.raises(subprocess.CalledProcessError):
        mock_sketch_obj.create_leaf_sketch()

def test_check_cardinality(mock_sketch_obj):
    """Test cardinality checking"""
    # Test when cardinality is not in cardkey
    mock_sketch_obj.sfp.full = "test_sketch"
    mock_sketch_obj.speciesinfo.cardkey = {}
    mock_sketch_obj._mock_sketch_check.return_value = True
    
    # Mock individual_card to set the cardkey value
    def mock_individual_card():
        mock_sketch_obj.speciesinfo.cardkey["test_sketch"] = 100.0
    mock_sketch_obj.individual_card = mock_individual_card
    
    card = mock_sketch_obj.check_cardinality()
    assert card == 100.0
    assert mock_sketch_obj.delta_pos == 100.0/21  # kval is 21
    
    # Test when cardinality is already in cardkey
    mock_sketch_obj.speciesinfo.cardkey["test_sketch"] = 200.0
    card = mock_sketch_obj.check_cardinality()
    assert card == 200.0
    assert mock_sketch_obj.delta_pos == 200.0/21 