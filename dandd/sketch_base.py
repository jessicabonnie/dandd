from __future__ import annotations
import os
import subprocess
import tempfile
from abc import ABC, abstractmethod
from dandd.species_specifics import SpeciesSpecifics
from dandd.sketch_filepath import SketchFilePath

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

    def create_sketch(self) -> None:
        '''Create a sketch or union of sketches based on the information used to initiate the sketch obj'''
        if not self.sketch_check():
            if self.ngen == 1:
                self.create_leaf_sketch()
            else:
                self.create_union_sketch()
            self.card = self.check_cardinality()
            self.delta_pos = self.card/self.kval

    def create_leaf_sketch(self) -> None:
        '''Create a sketch from a fasta file'''
        with tempfile.TemporaryDirectory() as tmpdir:
            self.cmd = self._leaf_command(tmpdir)
            if self.experiment['verbose']:
                print(self.cmd)
            try:
                subprocess.run(self.cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Command failed: {self.cmd}")
                raise e

    def create_union_sketch(self) -> None:
        '''Create a union of sketches'''
        self.cmd = self._union_command()
        if self.experiment['verbose']:
            print(self.cmd)
        try:
            subprocess.run(self.cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Command failed: {self.cmd}")
            raise e

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
            cmd = self.card_command()#[self.sfp.full])
            if not cmd:
                return
        stderr=None
        # if not self.experiment['verbose']:
        #     stderr=subprocess.DEVNULL
        if self.experiment['debug']:
            print(cmd)
        try:
            proc=subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, check=True,universal_newlines=True, stderr=stderr)
        except subprocess.CalledProcessError or RuntimeError:
            # warnings.warn(message=f"{self.sfp.full} cannot be created. Will attempt to remove and recreate component sketches.", category=RuntimeWarning)
            print(f"Recreating sketch {self.sfp.full}")
            self.create_sketch(just_do_it=True)
            # print("recreated sketch")
            proc=subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, universal_newlines=True, stderr=stderr)
        
        except subprocess.CalledProcessError as e:
            print(f"Command failed: {cmd}")
            raise e
        finally:
                self.parse_card(proc=proc)
                # NOTE: THIS IS NEW ... MAYBE IT BREAKS EVERYTHING SOON?
                self.check_cardinality()

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