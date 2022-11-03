from species_specifics import SpeciesSpecifics
import os
import hashlib
import warnings
import subprocess
import csv
DASHINGLOC="/home/jbonnie1/lib/dashing/dashing"

def canon_command(canon:bool):
        outstr=''
        if not canon:
            outstr=' --no-canon '
        return outstr
class SketchFilePath:
    '''An object to prepare sketch and union naming and directory location
    filenames: list of filenames that will be used to make the sketch
    kval: kvalue for the sketch
    speciesinfo: SpeciesSpecifics object for the species
    prefix: currently unused tag for filenames to differentiate between runs
    '''
    def __init__(self, filenames: list, kval: int, speciesinfo: SpeciesSpecifics, prefix=None):
        self.ffiles = filenames
        self.files= [os.path.basename(f) for f in self.ffiles]
        self.ngen = len(filenames)
        self.base = self.nameSketch(speciesinfo=speciesinfo, kval=kval)
        self.dir = os.path.join(speciesinfo.sketchdir,"k"+ str(kval), "ngen" + str(self.ngen))
        self.full = os.path.join(self.dir, self.base)
        self.registers = speciesinfo.registers
        os.makedirs(self.dir, exist_ok=True) 
    
    def assign_hash_string(self, filename: str, speciesinfo: SpeciesSpecifics, length: int):
        '''If the basename of filename already exists in hash/dict, look it up; otherwise create one and add to hash/dict using the full path as the key'''
        filename =os.path.basename(filename)
        if filename in speciesinfo.hashkey.keys():
            return speciesinfo.hashkey[filename]
        else:
            alphanum=hashlib.md5(filename.encode()).hexdigest()
            trunc=alphanum[:length]
            if trunc in speciesinfo.hashkey.values():
                warnings.warn("Hashvalue " + trunc + " has 2 keys!! " + filename + " will be assigned to a longer hash.")
                return self.assign_hash_string(filename, speciesinfo, length=length+1)
            else:
                speciesinfo.hashkey[filename] = trunc
                return trunc
    
    def nameSketch(self, speciesinfo: SpeciesSpecifics, kval: int):
        '''Determine what the name of a sketch is/will be'''
        if self.ngen > 1:
            
            filehashes = [self.assign_hash_string(filename=onefile, speciesinfo=speciesinfo, length=1) for onefile in self.files]
            filehashes.sort()
            outfile_prefix = speciesinfo.tag + "_" + "_".join(filehashes) + "_k" + str(kval) + "_r" + str(speciesinfo.registers)
            if len(outfile_prefix) >30:
                alphanum=hashlib.md5(outfile_prefix.encode()).hexdigest()
                outfile_prefix = alphanum
            outfile=outfile_prefix + ".hll"
        else:
            outfile=os.path.basename(self.files[0]) + ".w." + str(kval) + ".spacing." + str(speciesinfo.registers) + ".hll"
        return outfile    


class SketchObj:
    ''' A sketchobject.
        Inputs:
            kval: kmer length
            sfp: SketchFilePath object
            speciesinfo: SpeciesSpecifics object containing path information for the species/tag
            presketches: list of sketches to use in union sketching

        Properties:
            kval = the kvalue used to construct the sketch
            sketch = the location of the sketch file
            cmd = the command used to create the sketch
            card = the cardinality of the sketch
            delta_pos = possible delta value (to be compared to other ks of the same group of files)
    '''
    
    def __init__(self, kval, sfp, speciesinfo, presketches=None):
        self.kval = kval
        self.canon=canon_command(speciesinfo.canonicalize)
        self.sketch = None
        self.cmd = None
        self._sfp = sfp
        #self._registers = registers
        self._presketches = presketches
        #self._speciesinfo = speciesinfo
        self.create_sketch(sfp, speciesinfo)
        self.card = self.check_cardinality(speciesinfo)
        self.delta_pos = self.card/self.kval
    
    def __lt__(self, other):
        # lt = less than
        return self.delta_pos < other.delta_pos
    def __gt__(self, other):
        # gt = greater than
        return self.delta_pos > other.delta_pos
    def __repr__(self):
        ##TODO make this a __repr__ function instead
        return f"['sketch loc: {self.sketch}', k: {self.kval}, pos delta: {self.delta_pos}, cardinality: {self.card}, command: {self.cmd}  ]"
    
    
    def leaf_sketch(self, sfp, speciesinfo, debug=False):
        ''' If leaf sketch file exists, record the command that would have been used. If not run the command and store it.'''
        cmdlist = [DASHINGLOC, "sketch", "-k" + self.canon + str(self.kval),
                   "-S",str(sfp.registers),
                   "-p10","--prefix", str(sfp.dir),
                    sfp.ffiles[0]]
                #    os.path.join(speciesinfo.inputdir, sfp.files[0])]
        cmd = " ".join(cmdlist)
        if debug:
                print(cmd)
        if (not os.path.exists(sfp.full)) or os.stat(sfp.full).st_size == 0:
            print("The sketch file {0} either doesn't exist or is empty".format(sfp.full))
            subprocess.call(cmd, shell=True)
            self.cmd=cmd
            ##TODO Check if call raises error / returns 0
        else:
            self.cmd = cmd
        #print(self.cmd)
    #unionprefix <- file.path(sketchkndir,nameSketch(reorder[1:g], kval,registers=nregister))
     #command <- paste0("~/lib/dashing/dashing union -p ", parval," -z -o ", unionprefix, " ", alt_input1, " ", alt_input2)
        #print(command)
        #if (! file.exists(unionprefix) | file.size(unionprefix) == 0L){
        #  system(command, ignore.stdout = FALSE)
        #} 
        #sketch_call=subprocess.run([ "-p10","-o", str(sketchloc), left_sketch, right_sketch])
        #print(sketch_call)
    def union_sketch(self, sfp, debug=False):
        ''' If union sketch file exists, record the command that would have been used. If not run the command and store it.'''
        cmdlist = [DASHINGLOC, "union", "-p 10 ","-z -o", str(sfp.full)] + self._presketches
        cmd = " ".join(cmdlist)
        #print("SIZE of {0} is {1}".format(self._sfp.full, os.stat(sfp.full)))
        if (not os.path.exists(sfp.full)) or os.stat(sfp.full).st_size == 0:
            print("The sketch file {0} either doesn't exist or is empty".format(sfp.full))
            #self.cmd=subprocess.run(cmdlist)
            subprocess.call(cmd, shell=True)
            self.cmd=cmd
        else:
            self.cmd = cmd
            #" ".join(cmdlist)
        if debug:
            print(self.cmd)
    ##TODO: make leaf sketch and union sketch private
        
    def create_sketch(self, sfp: SketchFilePath, speciesinfo: SpeciesSpecifics):
        ''' If sketch file exists, assign path to self.sketch and return path. 
            If not create sketch, assign, and then return path.'''
        #self.sketch = self._sfp.full
        if sfp.ngen == 1:
            self.leaf_sketch(sfp, speciesinfo)
        elif sfp.ngen > 1:
            self.union_sketch(sfp)
        else:
            raise RuntimeError("For some reason you are trying to sketch an empty list of files. Don't do that.")
        self.sketch = sfp.full
        return self.sketch
    
    
    def individual_card(self, speciesinfo, debug=False):
        cmdlist = [DASHINGLOC,"card --presketched -p10"] +  [self.sketch]
        cmd = " ".join(cmdlist)
        if debug:
            print(cmd)
        #print(cmd)
        card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout.readlines()
        for card in csv.DictReader(card_lines, delimiter='\t'):
            speciesinfo.cardkey[card['#Path']] = card['Size (est.)']
            
    def check_cardinality(self, speciesinfo: SpeciesSpecifics):
        if self._sfp.full not in speciesinfo.cardkey.keys() or speciesinfo.cardkey[self._sfp.full] == 0:
            self.individual_card(speciesinfo)     
        return float(speciesinfo.cardkey[self.sketch])

        #else:
        #    self.individual_card(speciesinfo)
        #return float(speciesinfo.cardkey[self._sfp.full])
  