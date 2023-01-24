import pickle
import os
import re
import csv
from typing import List, Dict, Set, Tuple, NamedTuple

class SpeciesSpecifics:
    '''An object to store the specifics of a species file info'''
    def __init__(self, tag: str, genomedir: str, sketchdir: str, kstart: int, tool: str, flist_loc=None):
        self.tag=tag
        self.sketchdir=sketchdir
        #os.makedirs(self.sketchdir, exist_ok=True)
        self.species=self._resolve_species()
        self.fastahex=self._read_fastahex()
        self.cardkey=self._read_cardkey(tool=tool)
        self.inputdir=self._locate_input(genomedir)
        self.card0 = []
        self.kstart = kstart
        self.orderings = None
        self.flist_loc=flist_loc
        self.sketchinfo=self._read_sketchinfo()
    
    def read_pickle(self, filepath):
        if os.path.exists(filepath):
            contents=pickle.load(open(filepath, "rb", -1))
        else:
            contents=dict()
        return contents

    def _fastahex_loc(self):
        return os.path.join(self.sketchdir,'dandd_fastahex.pickle')
    def _sketchinfo_loc(self):
        # return os.path.join(self.sketchdir,self.tag+'_sketchinfo.pickle')
        return os.path.join(self.sketchdir,'dandd_sketchinfo.pickle')
    def _read_fastahex(self):
        '''Recover species specific fasta to hexidecimal dictionary from pickle file'''
        # usual=os.path.join(self.sketchdir,self.tag+'_fastahex.pickle')
        # usual=self._fastahex_loc()
        return self.read_pickle(self._fastahex_loc())
        # if os.path.exists(usual):
        #     fastahex=pickle.load(open(usual, "rb", -1))
        # else:
        #     fastahex=dict()
        # return fastahex
    def _read_sketchinfo(self) -> Dict[str,Dict]:
        '''Recover sketch name mappings from sketch directory file'''
        # usual=os.path.join(self.sketchdir,self.tag+'_fastahex.pickle')
        # usual=self._fastahex_loc()
        return self.read_pickle(self._sketchinfo_loc())
    ##TODO: Consider creating self.keys dict object to track the keys and writing a single save function
    # def save_key(self,obj,name):
    #     pass
    def _save_fastahex(self):
        '''Store/Update/Overwrite species specific hashkey to pickle'''
        # usual=os.path.join(self.sketchdir,'dandd_fastahex.pickle')
        # usual=os.path.join(self.sketchdir,self.tag+'_fastahex.pickle')
        with open(self._fastahex_loc(),"wb") as f:
            pickle.dump(file=f, obj=self.fastahex)
    
    def _save_sketchinfo(self):
        '''Store/Update/Overwrite sketchinfo lookup to pickle'''
        with open(self._sketchinfo_loc(),"wb") as f:
            pickle.dump(file=f, obj=self.sketchinfo)

    # def _write_experiment_key(self, bases:list, filepath:str):
    #     expdict = self.sketchinfo.fromkeys(bases)
    #     keys=self.sketchinfo[bases[0]].keys()
    #     with open(filepath, "w") as writer:
    #         dict_writer = csv.DictWriter(writer, fieldnames=keys)
    #         dict_writer.writeheader()
    #         dict_writer.writerows(expdict)
        # writer.close()
        
    def save_references(self, fast=True):
        if not fast:
            self._save_fastahex()
            self._save_sketchinfo()
        

    
    def _read_cardkey(self, tool):
        '''Recover key of previously calculated cardinalities from pickle file'''
        cardpath=os.path.join(self.sketchdir, f'{self.tag}_{tool}_cardinalities.pickle')
        # if os.path.exists(cardpath):
        #     cardkey=pickle.load(open(cardpath, "rb", -1))
        # else:
        #     cardkey=dict()
        return self.read_pickle(cardpath)
    
    def save_cardkey(self, tool: str, fast=False):
        '''Store cardinalities in species specific pickle'''
        if not fast:
            cardpath=os.path.join(self.sketchdir, f'{self.tag}_{tool}_cardinalities.pickle')
            with open(cardpath,"wb") as f:
                pickle.dump(file=f, obj=self.cardkey)
            
    def _resolve_species(self):
        '''Parse the name of the species from the data description tag'''
        for title in ['HVSVC2','ecoli','salmonella','human']:
            if title in self.tag:
                return title
        # os.warnings("FOR SOME REASON I DON'T RECOGNIZE THAT SPECIES TAG!!!")
        return self.tag
    
    def _locate_input(self, genomedir: str):
        '''Determine if there is a subdirectory structure in the input directory (current applies to HVSVC2 inputs)'''
        
        #grandfather in some old code
        if os.path.exists(os.path.join(genomedir, self.species)):
            if self.species == 'HVSVC2':
                return os.path.join(genomedir, self.species, self.tag.replace('HVSVC2','consensus'))
            return(os.path.join(genomedir, self.tag))
        else:
            return os.path.join(genomedir)
        
    def check_cardinality(self, fullpath):
        '''NOT IN USE. Check if cardinality has already been calculated for a sketch. If so, return it, if not add it to a list of cardinalities to be calculated and then return 0'''
        if fullpath in self.cardkey.keys() and self.cardkey[fullpath] != 0:
            return float(self.cardkey[fullpath])
        elif fullpath not in self.card0:
            self.card0.append(fullpath)
        return 0

    def retrieve_fasta_files(self, full=True)->list:
        '''return a list of all fasta files in a directory accounting for all the possible extensions'''
        reg_compile = re.compile(self.inputdir + "/*\.(fa.gz|fasta.gz|fna.gz|fasta|fa)")
        fastas = [fasta for fasta in os.listdir(self.inputdir) if reg_compile]
        if full:
            fastas=[os.path.join(self.inputdir,fasta) for fasta in fastas]
        return fastas