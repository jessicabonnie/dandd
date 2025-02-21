import pytest # type: ignore
import os
from dandd.utils import insert_pre_ext, write_listdict_to_csv, blake2b, canon_command, read_pickle_dict

def test_insert_pre_ext():
    assert insert_pre_ext("test.txt", "suffix") == "test.suffix.txt"
    assert insert_pre_ext("path/to/file.csv", "processed") == "path/to/file.processed.csv"
    assert insert_pre_ext("multiple.dots.file.txt", "new") == "multiple.dots.file.new.txt"

def test_write_listdict_to_csv(tmp_path):
    test_data = [
        {"name": "item1", "value": 10, "files": "file1.txt,file2.txt"},
        {"name": "item2", "value": 20, "files": "file3.txt"}
    ]
    
    output_file = os.path.join(tmp_path, "test_output.csv")
    write_listdict_to_csv(output_file, test_data)
    
    # Verify the file exists and contains expected headers
    assert os.path.exists(output_file)
    with open(output_file) as f:
        content = f.readlines()
        headers = content[0].strip().split(',')
        # Verify all expected headers are present
        assert all(h in headers for h in ['name', 'value', 'files'])
        # Verify 'files' is the last column (as per the function's design)
        assert headers[-1] == 'files'
        # Verify data is present
        assert 'item1' in content[1]
        assert 'file1.txt,file2.txt' in content[1]

def test_blake2b(tmp_path):
    # Create a test file with known content
    test_file = tmp_path / "test_file.txt"
    test_content = b"Hello, World!"
    test_file.write_bytes(test_content)
    
    # Calculate hash
    result = blake2b(str(test_file))
    
    # Verify it's a valid blake2b hash (128 hex characters)
    assert len(result) == 128
    # Verify it's consistent
    assert blake2b(str(test_file)) == result

def test_canon_command():
    # Test dashing tool
    assert canon_command(True, 'dashing') == ''
    assert canon_command(False, 'dashing') == '--no-canon'
    
    # Test kmc tool
    assert canon_command(True, 'kmc') == ''
    assert canon_command(False, 'kmc') == '-b'
    
    # Test default tool (dashing)
    assert canon_command(False) == '--no-canon' 

def test_read_pickle_dict(tmp_path):
    """Test pickle reading functionality"""
    # Test with non-existent file
    nonexistent = tmp_path / "nonexistent.pickle"
    result = read_pickle_dict(str(nonexistent))
    assert isinstance(result, dict)
    assert len(result) == 0

    # Test with valid pickle file
    test_data = {"key": "value"}
    test_pickle = tmp_path / "test.pickle"
    with open(test_pickle, "wb") as f:
        pickle.dump(test_data, f)
    
    result = temp_species.read_pickle(str(test_pickle))
    assert result == test_data