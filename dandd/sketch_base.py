from __future__ import annotations
import os
import subprocess
import tempfile
from abc import ABC, abstractmethod
from dandd.species_specifics import SpeciesSpecifics
from dandd.sketch_filepath import SketchFilePath
import shutil
DASHINGLOC="/home/jbonnie1/lib/dashing_working/dashing" 

class SketchObj(ABC):
    ''' Abstract base class for sketch objects '''
    
    def __init__(self, kval: int, sfp: SketchFilePath, speciesinfo: SpeciesSpecifics, experiment: dict, presketches=[]):
        self.kval = kval
        self.sketch = None
        self.cmd = None
        self.sfp = sfp
        self.delta_pos = 0
        self.card = 0
        self.speciesinfo = speciesinfo
        experiment['baseset'].add(sfp.base)
        self.experiment = experiment
        self._presketches = presketches
        if self.kval > 0:
            self.create_sketch()
            self.card = self.check_cardinality()
            self.delta_pos = self.card/self.kval

    def create_sketch(self, just_do_it=False):
        ''' If sketch file exists, assign path to self.sketch and return path. 
            If not create sketch, assign, and then return path.'''
        if self.sfp.ngen == 1:
            self.create_leaf_sketch(just_do_it=just_do_it)
        elif self.sfp.ngen > 1:
            self.create_union_sketch(just_do_it=just_do_it)
        else:
            raise RuntimeError("For some reason you are trying to sketch an empty list of files. Don't do that.")
        
        self.sketch = self.sfp.full
        return self.sketch

    def create_leaf_sketch(self, just_do_it=False):
        '''If leaf sketch file exists, record the command that would have been used. If not run the command and store it.'''
        tmpdir = tempfile.mkdtemp()
        cmd = self._leaf_command(tmpdir=tmpdir)
        stdout = None
        if not self.experiment['verbose']:
            stdout = subprocess.DEVNULL
        if self.experiment['debug']:
            print(cmd)
        try:
            if just_do_it or not self.sketch_check():
                if not (self.experiment["lowmem"] and self.check_cardinality() > 0):
                    if self.experiment["verbose"]:
                        print("Running Leaf Command: " + cmd)
                    try:
                        subprocess.run(cmd, shell=True, stdout=stdout, check=True)
                    except subprocess.CalledProcessError as e:
                        if e.returncode == 127:  # Command not found
                            tool = 'dashing' if 'dashing' in cmd else 'kmc'
                            raise RuntimeError(f"Could not find {tool} command. Please ensure it is installed and in your PATH.") from e
                        raise  # Re-raise other errors
                    self.cmd = cmd
            else:
                self.cmd = cmd
        finally:
            shutil.rmtree(tmpdir)

    def create_union_sketch(self, just_do_it=False):
        ''' If union sketch file exists, record the command that would have been used. If not run the command and store it.'''
        cmd = self._union_command()
        stdout=None
        stdout=subprocess.DEVNULL
        if self.experiment['verbose']:
            stdout=subprocess.PIPE
        if just_do_it:
            subprocess.call(cmd, shell=True, stdout=stdout)
            self.cmd = cmd
        elif not self.sketch_check():
            if self.experiment["lowmem"] and self.check_cardinality() > 0:
                self.cmd=cmd
            else:
                if self.experiment["verbose"]:
                    print("Running Union Command: " + cmd)
                subprocess.call(cmd, shell=True, stdout=stdout)
                self.cmd = cmd
        
        # self.cmd = cmd
        if self.experiment['debug']:
            print(self.cmd)

    def check_cardinality(self) -> float:
        '''Check whether the cardinality of sketch/db is stored in the cardkey, if not run a card command for the sketch and store it.'''
        if self.sfp.full not in self.speciesinfo.cardkey.keys():
            if self.experiment["lowmem"]:
                return 0

        if (self.sfp.full not in self.speciesinfo.cardkey.keys() or 
            float(self.speciesinfo.cardkey[self.sfp.full] == 0)) or self.speciesinfo.cardkey[self.sfp.full] is None:
            if not self.sketch_check():
                return 0
            self.individual_card()

        self.card = float(self.speciesinfo.cardkey[self.sfp.full])
        self.delta_pos = self.card/int(self.kval)
        return float(self.card)

    def individual_card(self, cmd=None) -> None:
        '''Run cardinality for an individual sketch or database. Add it to a dictionary {path:value}'''
        if self.kval == 0:
            return
        if not cmd:
            cmd = self.card_command([self.sfp.full])
            if not cmd:
                return
        stderr = None
        if self.experiment['debug']:
            print(cmd)
        try:
            proc = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, check=True, universal_newlines=True, stderr=stderr)
            self.parse_card(proc=proc)
            self.check_cardinality()
        except subprocess.CalledProcessError as e:
            if not self.experiment['mock_run']:  # Don't recreate if we're in test mode
                print(f"Recreating sketch {self.sfp.full}")
                self.create_sketch(just_do_it=True)
                proc = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, check=True, universal_newlines=True, stderr=stderr)
                self.parse_card(proc=proc)
                # NOTE: THIS IS NEW ... MAYBE IT BREAKS EVERYTHING SOON?
                self.check_cardinality()
            else:
                raise e

    def check_cardinality(self) -> float:
        '''Check whether the cardinality of sketch/db is stored in the cardkey, if not run a card command for the sketch and store it. '''
        if self.sfp.full not in self.speciesinfo.cardkey.keys():
            if self.experiment["lowmem"]:
                return 0
            
        # If the full path is not in the list of keys in the cardinality dictionary or the stored cardinality is 0, we will need to check if there is a sketch
        if (self.sfp.full not in self.speciesinfo.cardkey.keys() or float(self.speciesinfo.cardkey[self.sfp.full] == 0)) or self.speciesinfo.cardkey[self.sfp.full] is None:
            if not self.sketch_check():
                return 0
            self.individual_card()
        # Major error previously due to indentation
        self.card = float(self.speciesinfo.cardkey[self.sfp.full])
        self.delta_pos = self.card/int(self.kval)
        # if self.experiment["lowmem"] and self.sfp.ngen > 1:
        #     self.remove_sketch()
        return float(self.card)
    
    @abstractmethod
    def sketch_check(self, path=None) -> bool:
        """Check if sketch exists and is valid"""
        pass

    @abstractmethod
    def _leaf_command(self, tmpdir) -> str:
        """Generate command to create leaf sketch"""
        pass

    @abstractmethod
    def remove_sketch(self, delete_me: str = None):
        """Remove sketch file(s)"""
        pass

    @abstractmethod
    def card_command(self, sketch_paths=[]) -> str:
        """Generate command to get cardinality"""
        pass

    @abstractmethod
    def parse_card(self, proc):
        """Parse cardinality command output"""
        pass

    @abstractmethod
    def _union_command(self) -> str:
        """Generate command to create union of sketches"""
        pass

    @property
    def ngen(self) -> int:
        """Number of input files"""
        return self.sfp.ngen 