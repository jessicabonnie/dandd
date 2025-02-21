#!/usr/bin/env python3
import os
import pickle
import sqlite3
import argparse
from pathlib import Path
import sys

## NOTE CURRENTLY NOT USED or TESTED

def migrate_cardinalities(pickle_dir: str, db_path: str, verbose: bool = False) -> None:
    """Migrate cardinality data from pickle files to SQLite database"""
    
    # Initialize database and create table
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS cardinalities (
                sketch_path TEXT PRIMARY KEY,
                tool TEXT NOT NULL,
                tag TEXT NOT NULL,
                cardinality REAL NOT NULL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        conn.commit()

    # Find all cardinality pickle files
    pickle_files = []
    for file in Path(pickle_dir).rglob("*_cardinalities.pickle"):
        if file.is_file():
            pickle_files.append(file)

    if verbose:
        print(f"Found {len(pickle_files)} cardinality pickle files")

    # Process each pickle file
    for pickle_file in pickle_files:
        try:
            # Extract tag and tool from filename
            filename = pickle_file.stem  # Get filename without extension
            parts = filename.split('_')
            if len(parts) >= 3 and parts[-1] == 'cardinalities':
                tag = parts[0]
                tool = parts[1]
            else:
                print(f"Warning: Unexpected filename format: {pickle_file}", file=sys.stderr)
                continue

            # Load pickle data
            with open(pickle_file, 'rb') as f:
                try:
                    cardkey = pickle.load(f)
                except (pickle.UnpicklingError, EOFError):
                    # Try backup file if it exists
                    backup_file = pickle_file.with_suffix('.pickle.bkp')
                    if backup_file.exists():
                        with open(backup_file, 'rb') as bf:
                            try:
                                cardkey = pickle.load(bf)
                            except (pickle.UnpicklingError, EOFError):
                                print(f"Error: Both {pickle_file} and backup are corrupted", file=sys.stderr)
                                continue
                    else:
                        print(f"Error: {pickle_file} is corrupted and no backup exists", file=sys.stderr)
                        continue

            # Insert data into SQLite
            with sqlite3.connect(db_path) as conn:
                cursor = conn.cursor()
                cursor.execute('BEGIN TRANSACTION')
                try:
                    for sketch_path, cardinality in cardkey.items():
                        cursor.execute('''
                            INSERT OR REPLACE INTO cardinalities 
                            (sketch_path, tool, tag, cardinality)
                            VALUES (?, ?, ?, ?)
                        ''', (sketch_path, tool, tag, cardinality))
                    conn.commit()
                    if verbose:
                        print(f"Migrated {len(cardkey)} entries from {pickle_file}")
                except Exception as e:
                    conn.rollback()
                    print(f"Error migrating {pickle_file}: {str(e)}", file=sys.stderr)
                    continue

        except Exception as e:
            print(f"Error processing {pickle_file}: {str(e)}", file=sys.stderr)
            continue

def main():
    parser = argparse.ArgumentParser(description='Migrate DandD cardinality data from pickle files to SQLite')
    parser.add_argument('pickle_dir', help='Directory containing pickle files')
    parser.add_argument('--db-path', default='dandd.db', help='Path to SQLite database (default: dandd.db)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose output')
    args = parser.parse_args()

    if not os.path.isdir(args.pickle_dir):
        print(f"Error: {args.pickle_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    try:
        migrate_cardinalities(args.pickle_dir, args.db_path, args.verbose)
        print("Migration completed successfully")
    except Exception as e:
        print(f"Migration failed: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()