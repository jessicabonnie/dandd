import pickle
import os
import re
from typing import Dict
import shutil
from subprocess import CalledProcessError

class SpeciesSpecifics:
    '''An object to store the specifics of a species file info'''
    def __init__(self, tag: str, genomedir: str, sketchdir: str, kstart: int, tool: str, flist_loc=None):
        self.tag=tag
        self.sketchdir=sketchdir
        # self.species=self._resolve_species()
        self.fastahex=self._read_fastahex()
        self.cardkey=self._read_cardkey(tool=tool)
        self.inputdir=genomedir
        self.card0 = []
        self.kstart = kstart
        self.orderings = None
        self.flist_loc=flist_loc
        self.sketchinfo=self._read_sketchinfo()
    
    def read_pickle(self, filepath) -> Dict:
        '''Read a pickle into a dictionary if the filepath exists, otherwise return an empty dictionary'''
        if os.path.exists(filepath):
            try:
                contents=pickle.load(open(filepath, "rb", -1))
            except CalledProcessError:
                contents=pickle.load(open(filepath+'.bkp', "rb", -1))
            # finally:
            #     print(f"{filepath} and {filepath}.bkp are both corrupted. They will be overwritten.")
            #     contents=dict()
        else:
            contents=dict()
        return contents

    def _fastahex_loc(self)-> str:
        return os.path.join(self.sketchdir,'dandd_fastahex.pickle')
    def _sketchinfo_loc(self)-> str:
        return os.path.join(self.sketchdir,'dandd_sketchinfo.pickle')
    def _read_fastahex(self):
        '''Recover species specific fasta to hexidecimal dictionary from pickle file'''
        return self.read_pickle(self._fastahex_loc())

    def _read_sketchinfo(self) -> Dict[str,Dict]:
        '''Recover sketch name mappings from sketch directory file'''
        return self.read_pickle(self._sketchinfo_loc())

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

    
    def _read_cardkey(self, tool) -> Dict[str, int]:
        '''Recover key of previously calculated cardinalities from pickle file'''
        cardpath=os.path.join(self.sketchdir, f'{self.tag}_{tool}_cardinalities.pickle')
        return self.read_pickle(cardpath)
    
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