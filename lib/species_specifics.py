import pickle
import os

class SpeciesSpecifics:
    '''An object to store the specifics of a species file info'''
    def __init__(self, tag: str, genomedir: str, sketchdir: str, kstart: int, registers=20):
        self.tag=tag
        self.sketchdir=os.path.join(sketchdir, tag)
        os.makedirs(self.sketchdir, exist_ok=True)
        self.species=self._resolve_species()
        self.hashkey=self._read_hashkey()
        self.cardkey=self._read_cardkey()
        self.inputdir=self._locate_input(genomedir)
        self.card0 = []
        self.kstart = kstart
        self.orderings = None
        self.registers=registers
        
    def _read_hashkey(self):
        '''Recover species specific hashkey from pickle file'''
        usual=os.path.join(self.sketchdir,self.tag+'_hashkey.pickle')
        if os.path.exists(usual):
            hashkey=pickle.load(open(usual, "rb", -1))
        else:
            hashkey=dict()
        return hashkey

    def save_hashkey(self):
        '''Store/Update/Overwrite species specific hashkey to pickle'''
        usual=os.path.join(self.sketchdir,self.tag+'_hashkey.pickle')
        with open(usual,"wb") as f:
            pickle.dump(file=f, obj=self.hashkey)
    
    def _read_cardkey(self):
        '''Recover key of previously calculated cardinalities from pickle file'''
        usual=os.path.join(self.sketchdir, self.tag+'_cardinalities.pickle')
        if os.path.exists(usual):
            cardkey=pickle.load(open(usual, "rb", -1))
        else:
            cardkey=dict()
        return cardkey
    
    def save_cardkey(self):
        '''Store cardinalities in species specific pickle'''
        usual=os.path.join(self.sketchdir, self.tag+'_cardinalities.pickle')
        with open(usual,"wb") as f:
            pickle.dump(file=f, obj=self.cardkey)
            
    def _resolve_species(self):
        '''Parse the name of the species from the data description tag'''
        for title in ['HVSVC2','ecoli','salmonella','human']:
            if title in self.tag:
                return title
        os.warnings("FOR SOME REASON I DON'T RECOGNIZE THAT SPECIES TAG!!!")
        return self.tag
    
    def _locate_input(self, genomedir: str):
        '''Determine if there is a subdirectory structure in the input directory (current applies to HVSVC2 inputs)'''
        if self.species == 'HVSVC2':
            return os.path.join(genomedir, self.species, self.tag.replace('HVSVC2','consensus'))
        else:
            return os.path.join(genomedir, self.tag)
        
    def check_cardinality(self, fullpath):
        '''Check if cardinality has already been calculated for a sketch. If so, return it, if not add it to a list of cardinalities to be calculated and then return 0'''
        if fullpath in self.cardkey.keys() and self.cardkey[fullpath] != 0:
            return float(self.cardkey[fullpath])
        elif fullpath not in self.card0:
            self.card0.append(fullpath)
        return 0