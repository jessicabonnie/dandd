import pytest
import os
import subprocess
from unittest.mock import Mock, patch
from dandd.species_specifics import SpeciesSpecifics
from dandd.utils import blake2b

@pytest.fixture
def mock_run():
    """Mock subprocess.run"""
    with patch('subprocess.run') as mock:
        mock.return_value = Mock(returncode=0)
        yield mock

@pytest.fixture
def test_experiment():
    """Basic experiment configuration"""
    return {
        'registers': 16,
        'canonicalize': True,
        'tool': 'dashing',
        'safety': True,
        'verbose': False,
        'debug': False,
        'lowmem': False,
        'nthreads': 1,
        'baseset': set()
    }

@pytest.fixture
def test_data_dir():
    """Path to test data directory"""
    return os.path.join(os.path.dirname(__file__), 'data')

@pytest.fixture
def test_files(test_data_dir):
    """Get test files with proper paths"""
    return [
        os.path.join(test_data_dir, "NC_009057.fasta"),
        os.path.join(test_data_dir, "NC_009058.fasta.gz")
    ]

@pytest.fixture
def temp_species(tmp_path, test_data_dir):
    """Create a temporary SpeciesSpecifics instance with test directories"""
    genomedir = tmp_path / "genomes"
    sketchdir = tmp_path / "sketches"
    genomedir.mkdir()
    sketchdir.mkdir()
    
    # Copy test files to genomedir
    import shutil
    for fname in ["NC_009057.fasta", "NC_009058.fasta.gz"]:
        shutil.copy2(
            os.path.join(test_data_dir, fname),
            os.path.join(genomedir, fname)
        )
    
    species = SpeciesSpecifics(
        tag="test_species",
        genomedir=str(genomedir),
        sketchdir=str(sketchdir),
        kstart=21,
        tool="dashing"
    )
    
    # Initialize fastahex with actual blake2b hashes of test files
    species.fastahex = {}
    for fname in ["NC_009057.fasta", "NC_009058.fasta.gz"]:
        species.fastahex[fname] = blake2b(os.path.join(genomedir, fname))
    
    # Initialize cardkey
    species.cardkey = {}
    
    # Create sketch directory
    os.makedirs(species.sketchdir, exist_ok=True)
    
    return species 