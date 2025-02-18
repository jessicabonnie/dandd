from __future__ import annotations
import os
from dandd.species_specifics import SpeciesSpecifics
from dandd.utils import blake2b

class SketchFilePath:
    '''An object to prepare sketch and union naming and directory location
    filenames: list of filenames that will be used to make the sketch
    kval: kmer length
    speciesinfo: SpeciesSpecifics object for the species
    prefix: currently unused tag for filenames to differentiate between runs
    '''
    def __init__(self, filenames: list, kval: int, speciesinfo: SpeciesSpecifics, experiment:dict, prefix=None ):
        self.ffiles = filenames
        self.files= [os.path.basename(f) for f in self.ffiles]
        self.files.sort()
        self.ngen = len(filenames)
        # If kval is 0 then we know we are in a ksweep situation
        kval_str=str(kval)
        if kval == 0:
            kval_str="{}"
        self.dir = os.path.join(speciesinfo.sketchdir, "ngen" + str(self.ngen),"k"+ kval_str)
        self.base =self._assign_base(speciesinfo=speciesinfo, kval=kval, registers=experiment['registers'], canonicalize=experiment['canonicalize'], tool=experiment['tool'], safety=experiment['safety'])
        self.relative = os.path.join("ngen" + str(self.ngen),"k"+ str(kval),self.base)+ self._get_ext(experiment['tool'])
        self.full = os.path.join(self.dir, self.base)+ self._get_ext(experiment['tool'])
        if kval != 0:
            os.makedirs(self.dir, exist_ok=True)

    def _get_ext(self,tool)->str:
        if tool == 'dashing':
            ext='.hll'
        elif tool == 'kmc':
            ext=''
        else:
            raise ValueError("is there another option for tool other than kmc or dashing?")
        return ext
        
    def _hashsum(self, speciesinfo:SpeciesSpecifics):
        '''Calculate the blake2b hexsum of an individual fastas or sum the hexsums of component fastas to create hexidecimal identifiers for combinations of fastas'''
        if self.ngen == 1:
            # Return just the hash string without 0x prefix for single files
            output = blake2b(self.ffiles[0])
        else:
            sum = int("0",16)
            for fasta in self.files:
                if fasta not in speciesinfo.fastahex:
                    speciesinfo.fastahex[fasta] = blake2b(os.path.join(speciesinfo.inputdir, fasta))
                # Remove 0x prefix if present when converting to int
                hash_val = speciesinfo.fastahex[fasta]
                if hash_val.startswith('0x'):
                    hash_val = hash_val[2:]
                sum += int(hash_val, 16)
            output = hex(sum)[2:]  # Remove 0x prefix from final sum
        return output
        
    def _assign_base(self, speciesinfo:SpeciesSpecifics, kval:int, registers:int, canonicalize:bool, tool:str, safety=False) -> str:
        '''determine the base file name for the sketch using the properties that will be used to generate it'''
        kval_str=str(kval)
        if kval == 0:
            kval_str="{}"
        fnames_key=''.join(self.files)
        # if the key (made by joining the ingredient filenames) isn't already in the fastahex dictionary mapping the combination of those files to a hexsum, calculate that hexsum and add it to the fastahex key
        if fnames_key not in speciesinfo.fastahex.keys():
            speciesinfo.fastahex[fnames_key]= self._hashsum(speciesinfo)   
            stored_val=speciesinfo.fastahex[fnames_key]
        # if the key is there, calculate what we expect the hashsum value to be based on the hexes of the components -- this is just to confirm that nothing has gotten confused somehow
        else:
            stored_val=speciesinfo.fastahex[fnames_key]
            if safety:
                checkval=self._hashsum(speciesinfo)
                if checkval != stored_val:
                    raise RuntimeError(f"Checksum does not match stored value for {fnames_key}: {checkval}, {stored_val}")
        # if the sketch is of a single input file, dashing will insist on naming it something specific, so we will use that base as a name for both dashing and kmc to make life easier
        if self.ngen == 1:
            if tool == 'dashing':
                sketchbase=self.files[0] + ".w." + kval_str + ".spacing." + str(registers)
            else:
                sketchbase=self.files[0] + "_k" + kval_str
                if not canonicalize:
                    sketchbase=sketchbase+'nc'
        # if sketch is of a combination of sketches, add a tag that will differentiate it from other combinations of the same sketches composed using different register counts or kvalues. Also guard against the low probability chance that there are overlapping hexsums of different numbers of fasta inputs
        else:
            suffix = str(registers) + "n" + str(self.ngen) + "k" + kval_str
            if not canonicalize:
                suffix=suffix+'nc'
            sketchbase = stored_val[:15] + "_" + suffix
        # store information relating to this basename to be given to user later as table or obj
        info = {"sketchbase": sketchbase, "files": self.files, "ngen": self.ngen, "kval": kval, "registers": registers }
        if sketchbase not in speciesinfo.sketchinfo.keys():
            speciesinfo.sketchinfo[sketchbase] = info
        elif safety:
            stored_info=speciesinfo.sketchinfo[sketchbase]
            # check to make sure all sketchinfo values match what is stored
            for key in stored_info.keys():
              if stored_info[key] != info[key]:
                  raise RuntimeError(f"Duplicate keys but not duplicate values: {sketchbase}: (1) {stored_info}, (2) {info}")
        return sketchbase