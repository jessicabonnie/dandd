import pytest # type: ignore
import os
from dandd.sketch_filepath import SketchFilePath
from dandd.species_specifics import SpeciesSpecifics
from dandd.utils import blake2b

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
def temp_species(tmp_path):
    """Create a temporary SpeciesSpecifics instance with test directories"""
    genomedir = tmp_path / "genomes"
    sketchdir = tmp_path / "sketches"
    genomedir.mkdir()
    sketchdir.mkdir()
    
    # Create test fasta files
    (genomedir / "test1.fa").write_text(">seq1\nACGT\n")
    (genomedir / "test2.fasta.gz").write_text(">seq2\nGCTA\n")
    
    species = SpeciesSpecifics(
        tag="test_species",
        genomedir=str(genomedir),
        sketchdir=str(sketchdir),
        kstart=21,
        tool="dashing"
    )
    
    # Initialize fastahex with actual blake2b hashes
    species.fastahex = {
        "test1.fa": blake2b(str(genomedir / "test1.fa")),
        "test2.fasta.gz": blake2b(str(genomedir / "test2.fasta.gz"))
    }
    
    return species

@pytest.fixture
def test_files(temp_species):
    """Get test files with proper paths"""
    return [
        os.path.join(temp_species.inputdir, "test1.fa"),
        os.path.join(temp_species.inputdir, "test2.fasta.gz")
    ]

def test_init(temp_species, test_experiment, test_files):
    """Test initialization of SketchFilePath"""
    sfp = SketchFilePath(
        filenames=test_files,
        kval=21,
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    
    # Debug print
    print(f"base: {sfp.base}")
    print(f"files: {sfp.files}")
    print(f"fastahex: {temp_species.fastahex}")
    
    assert sfp.ngen == 2
    assert sfp.files == ["test1.fa", "test2.fasta.gz"]  # Sorted order
    assert sfp.base is not None  # Make sure base is assigned
    assert sfp.dir.endswith(os.path.join("ngen2", "k21"))
    assert sfp.full.endswith(".hll")  # For dashing tool

def test_get_ext(temp_species, test_experiment, test_files):
    """Test file extension generation"""
    sfp = SketchFilePath(
        filenames=[test_files[0]],  # Use single file to simplify test
        kval=21,
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    assert sfp._get_ext('dashing') == '.hll'
    assert sfp._get_ext('kmc') == ''
    with pytest.raises(ValueError):
        sfp._get_ext('invalid_tool')

def test_hashsum(temp_species, test_experiment, test_files):
    """Test hash sum calculation"""
    # Test single file
    sfp_single = SketchFilePath(
        filenames=[test_files[0]],  # Single file
        kval=21,
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    
    single_hash = sfp_single._hashsum(temp_species)
    assert isinstance(single_hash, str)
    # We don't expect 0x prefix anymore based on implementation
    assert len(single_hash) == 128  # blake2b produces 64-byte (128 hex char) hashes
    
    # Test multiple files
    sfp_multi = SketchFilePath(
        filenames=test_files,  # Multiple files
        kval=21,
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    
    # Store first hash for later comparison
    first_hash = blake2b(test_files[0])
    temp_species.fastahex[os.path.basename(test_files[0])] = first_hash
    
    # Store second hash
    second_hash = blake2b(test_files[1])
    temp_species.fastahex[os.path.basename(test_files[1])] = second_hash
    
    multi_hash = sfp_multi._hashsum(temp_species)
    assert isinstance(multi_hash, str)
    # Verify we can convert both hashes to integers and sum them
    sum_of_hashes = hex(int(first_hash, 16) + int(second_hash, 16))[2:]  # Remove 0x prefix
    assert multi_hash == sum_of_hashes

def test_assign_base_single_file(temp_species, test_experiment, test_files):
    """Test base name assignment for single file"""
    sfp = SketchFilePath(
        filenames=[test_files[0]],
        kval=21,
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    
    base = sfp._assign_base(
        speciesinfo=temp_species,
        kval=21,
        registers=test_experiment['registers'],
        canonicalize=test_experiment['canonicalize'],
        tool=test_experiment['tool'],
        safety=test_experiment['safety']
    )
    
    assert base is not None
    assert "w.21.spacing.16" in base  # For dashing
    assert "test1.fa" in base

def test_assign_base_multiple_files(temp_species, test_experiment, test_files):
    """Test base name assignment for multiple files"""
    sfp = SketchFilePath(
        filenames=test_files,
        kval=21,
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    
    base = sfp._assign_base(
        speciesinfo=temp_species,
        kval=21,
        registers=test_experiment['registers'],
        canonicalize=test_experiment['canonicalize'],
        tool=test_experiment['tool'],
        safety=test_experiment['safety']
    )
    
    assert base is not None
    assert "n2k21" in base  # Multiple files indicator
    assert len(base) > 15  # Should include hash prefix

def test_ksweep_handling(temp_species, test_experiment, test_files):
    """Test handling of k-sweep mode (kval=0)"""
    sfp = SketchFilePath(
        filenames=test_files,
        kval=0,  # k-sweep mode
        speciesinfo=temp_species,
        experiment=test_experiment
    )
    
    assert sfp.base is not None
    assert "k{}" in sfp.dir
    assert not os.path.exists(sfp.dir)  # Directory shouldn't be created for k=0 