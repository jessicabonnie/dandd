import pytest # type: ignore
import os
import pickle
from dandd.species_specifics import SpeciesSpecifics
from dandd.utils import read_pickle_dict

@pytest.fixture
def temp_species(tmp_path):
    """Create a temporary SpeciesSpecifics instance with test directories"""
    genomedir = tmp_path / "genomes"
    sketchdir = tmp_path / "sketches"
    genomedir.mkdir()
    sketchdir.mkdir()
    
    # Create some test fasta files
    (genomedir / "test1.fa").write_text("test sequence 1")
    (genomedir / "test2.fasta.gz").write_text("test sequence 2")
    
    return SpeciesSpecifics(
        tag="test_species",
        genomedir=str(genomedir),
        sketchdir=str(sketchdir),
        kstart=21,
        tool="dashing"
    )

def test_init(temp_species):
    """Test initialization of SpeciesSpecifics"""
    assert temp_species.tag == "test_species"
    assert temp_species.kstart == 21
    assert isinstance(temp_species.fastahex, dict)
    assert isinstance(temp_species.cardkey, dict)
    assert isinstance(temp_species.sketchinfo, dict)

# def test_read_pickle(temp_species, tmp_path):
#     """Test pickle reading functionality"""
#     # Test with non-existent file
#     nonexistent = tmp_path / "nonexistent.pickle"
#     result = temp_species.read_pickle(str(nonexistent))
#     assert isinstance(result, dict)
#     assert len(result) == 0

#     # Test with valid pickle file
#     test_data = {"key": "value"}
#     test_pickle = tmp_path / "test.pickle"
#     with open(test_pickle, "wb") as f:
#         pickle.dump(test_data, f)
    
#     result = temp_species.read_pickle(str(test_pickle))
#     assert result == test_data

def test_retrieve_fasta_files(temp_species):
    """Test fasta file retrieval"""
    files = temp_species.retrieve_fasta_files(full=True)
    assert len(files) == 2
    assert any(f.endswith('.fa') for f in files)
    assert any(f.endswith('.fasta.gz') for f in files)
    
    # Test without full paths
    short_files = temp_species.retrieve_fasta_files(full=False)
    assert len(short_files) == 2
    assert all(not os.path.isabs(f) for f in short_files)

def test_save_and_read_references(temp_species):
    """Test saving and reading reference data"""
    # Add some test data
    temp_species.fastahex = {"test": "hex123"}
    temp_species.sketchinfo = {"sketch1": {"info": "data"}}
    
    # Save references
    temp_species.save_references()
    
    # Create new instance to read saved data
    new_species = SpeciesSpecifics(
        tag=temp_species.tag,
        genomedir=temp_species.inputdir,
        sketchdir=temp_species.sketchdir,
        kstart=temp_species.kstart,
        tool="dashing"
    )
    
    assert new_species.fastahex == temp_species.fastahex
    assert new_species.sketchinfo == temp_species.sketchinfo

def test_save_and_read_cardkey(temp_species):
    """Test saving and reading cardinality data"""
    # Add some test cardinality data
    temp_species.cardkey = {"sequence1": 100, "sequence2": 200}
    
    # Save cardkey
    temp_species.save_cardkey("dashing")
    
    # Create new instance to read saved data
    new_species = SpeciesSpecifics(
        tag=temp_species.tag,
        genomedir=temp_species.inputdir,
        sketchdir=temp_species.sketchdir,
        kstart=temp_species.kstart,
        tool="dashing"
    )
    
    assert new_species.cardkey == temp_species.cardkey
