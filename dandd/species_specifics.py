#! /usr/bin/env python3
"""
Module for managing species-specific data storage and retrieval.
Handles both SQLite database and pickle-based storage of:
- Fasta file hexsums
- Sketch cardinalities
- Sketch information
"""
from __future__ import annotations
import pickle
import os
import re
from typing import Dict, List, Optional
import shutil
from subprocess import CalledProcessError
from dandd.utils import read_pickle_dict
import sqlite3

class SpeciesSpecifics:
    """
    Manages species-specific data storage and retrieval.
    
    Attributes:
        tag (str): Identifier for the species/experiment
        sketchdir (str): Directory where sketches and metadata are stored
        fastahex (Dict[str, str]): Maps fasta paths to their hexsums
        db_path (str): Path to SQLite database
        cardkey (Dict[str, float]): Maps sketch paths to their cardinalities
        inputdir (str): Directory containing input fasta files
        card0 (List[str]): List of zero cardinality sketches
        kstart (int): Starting k-mer size
        orderings (Optional[List[List[str]]]): Ordering information for sketches
        flist_loc (Optional[str]): Location of fasta list file
        sketchinfo (Dict[str, Dict]): Mapping of sketch metadata
    """

    def __init__(self, tag: str, genomedir: str, sketchdir: str, kstart: int, tool: str, flist_loc: Optional[str] = None) -> None:
        """
        Initialize SpeciesSpecifics object.
        
        Args:
            tag: Identifier for the species/experiment
            genomedir: Directory containing input fasta files
            sketchdir: Directory where sketches and metadata are stored
            kstart: Starting k-mer size
            tool: Sketching tool to use ('dashing' or 'kmc')
            flist_loc: Optional path to file listing fasta files to use
        """
        self.tag=tag
        self.sketchdir=sketchdir
        # self.species=self._resolve_species()
        self.fastahex=self._read_fastahex()
        self.db_path = os.path.join(self.sketchdir, 'dandd.db')
        self._init_db()
        self.cardkey=self._read_cardkey(tool=tool)
        self.inputdir=genomedir
        self.card0 = []
        self.kstart = kstart
        self.orderings = None
        self.flist_loc=flist_loc
        self.sketchinfo=self._read_sketchinfo()
    

    def _fastahex_loc(self)-> str:
        return os.path.join(self.sketchdir,'dandd_fastahex.pickle')
    def _sketchinfo_loc(self)-> str:
        return os.path.join(self.sketchdir,'dandd_sketchinfo.pickle')
    def _read_fastahex(self):
        '''Recover species specific fasta to hexidecimal dictionary from pickle file'''
        return read_pickle_dict(self._fastahex_loc())

    def _read_sketchinfo(self) -> Dict[str,Dict]:
        '''Recover sketch name mappings from sketch directory file'''
        return read_pickle_dict(self._sketchinfo_loc())

    def update(self, tool) -> None:
        self.fastahex=self._read_fastahex()
        self.cardkey=self._read_cardkey(tool=tool)
        self.sketchinfo=self._read_sketchinfo()

    def _save_fastahex(self) -> None:
        '''Store/Update/Overwrite species specific hashkey to pickle'''
        fasta_hex_loc=self._fastahex_loc()
        with open(fasta_hex_loc+'.bkp',"wb") as f:
            pickle.dump(file=f, obj=self.fastahex)
        shutil.copy(fasta_hex_loc+'.bkp',fasta_hex_loc)
    
    def _save_sketchinfo(self) -> None:
        '''Store/Update/Overwrite sketchinfo lookup to pickle'''
        sketchinfo_loc=self._sketchinfo_loc()
        with open(sketchinfo_loc+'.bkp',"wb") as f:
            pickle.dump(file=f, obj=self.sketchinfo)
        shutil.copy(sketchinfo_loc+'.bkp', sketchinfo_loc)

    def save_references(self, fast=False) -> None:
        '''Save the fastahex and the sketchinfo objects to their default locations, overwriting what was there.'''
        if not fast:
            self._save_fastahex()
            self._save_sketchinfo()  

    
    def _read_cardkey(self, tool) -> Dict[str, float]:
        '''Recover key of previously calculated cardinalities from pickle file'''
        cardpath=os.path.join(self.sketchdir, f'{self.tag}_{tool}_cardinalities.pickle')
        return read_pickle_dict(cardpath)
    def _read_cardkey_db(self, tool: str) -> Dict[str, float]:
        """Read cardinalities from SQLite database"""
        cardkey = {}
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT sketch_path, cardinality 
                FROM cardinalities 
                WHERE tool = ? AND tag = ?
            ''', (tool, self.tag))
            for sketch_path, cardinality in cursor.fetchall():
                cardkey[sketch_path] = cardinality
        return cardkey

    def save_cardkey_db(self, tool: str, fast=False) -> None:
        """Store cardinalities in SQLite database"""
        if not fast:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                # Use a transaction for better performance with multiple inserts
                cursor.execute('BEGIN TRANSACTION')
                try:
                    for sketch_path, cardinality in self.cardkey.items():
                        cursor.execute('''
                            INSERT OR REPLACE INTO cardinalities 
                            (sketch_path, tool, tag, cardinality)
                            VALUES (?, ?, ?, ?)
                        ''', (sketch_path, tool, self.tag, cardinality))
                    conn.commit()
                except Exception as e:
                    conn.rollback()
                    raise e
    
    def save_cardkey(self, tool: str, fast=False) -> None:
        '''Store cardinalities in species specific pickle'''
        if not fast:
            cardpath=os.path.join(self.sketchdir, f'{self.tag}_{tool}_cardinalities.pickle')
            with open(cardpath+'.bkp',"wb") as f:
                pickle.dump(file=f, obj=self.cardkey)
            shutil.copyfile( cardpath+'.bkp', cardpath)

    def retrieve_fasta_files(self, full=True)->list:
        '''return a list of all fasta files in a directory accounting for all the possible extensions'''
        reg_compile = re.compile(self.inputdir + "/*\.(fa.gz|fasta.gz|fna.gz|fasta|fa)")
        fastas = [fasta for fasta in os.listdir(self.inputdir) if reg_compile]
        if full:
            fastas=[os.path.join(self.inputdir,fasta) for fasta in fastas]
        return fastas

    def _init_db(self):
        """Initialize SQLite database with necessary tables if they don't exist"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            # Cardinalities table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS cardinalities (
                    sketch_path TEXT PRIMARY KEY,
                    tool TEXT NOT NULL,
                    tag TEXT NOT NULL,
                    cardinality REAL NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            # Fastahex table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS fastahex (
                    fasta_path TEXT PRIMARY KEY,
                    hex_value TEXT NOT NULL,
                    tag TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            conn.commit()

    def save_fastahex_db(self, fast=False) -> None:
        """Store fasta hexsums in SQLite database"""
        if not fast:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute('BEGIN TRANSACTION')
                try:
                    for fasta_path, hex_value in self.fastahex.items():
                        cursor.execute('''
                            INSERT OR REPLACE INTO fastahex 
                            (fasta_path, hex_value, tag)
                            VALUES (?, ?, ?)
                        ''', (fasta_path, hex_value, self.tag))
                    conn.commit()
                except Exception as e:
                    conn.rollback()
                    raise e
    def _read_fastahex_db(self) -> Dict[str, str]:
        """Read fasta hexsums from SQLite database"""
        fastahex = {}
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT fasta_path, hex_value 
                FROM fastahex 
                WHERE tag = ?
            ''', (self.tag,))
            for fasta_path, hex_value in cursor.fetchall():
                fastahex[fasta_path] = hex_value
        return fastahex