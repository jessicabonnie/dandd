import os
import tempfile
import subprocess
import shutil
from typing import List
from dandd.species_specifics import SpeciesSpecifics
from dandd.sketch_dashing import DashSketchObj
from dandd.sketch_kmc import KMCSketchObj
from dandd.sketch_filepath import SketchFilePath

class DeltaNode:
    ''' A node in a Delta tree. 
        node_title = name of input file or composite of inputfiles
        children = the nodes that are this nodes children
        progeny = list of leaf nodes decended from the node
        ksweep = should ks 1-100 be explored even after delta is found?
        
        '''
    def __init__(self, node_title: str, children: list, speciesinfo: SpeciesSpecifics, experiment:dict, progeny:list=[]):
        self.node_title = node_title
        self.progeny = progeny
        self.experiment=experiment
        self.speciesinfo=speciesinfo
        self.children=children
        self.mink = 0
        self.maxk = 0
        if self.experiment["ksweep"] is not None:
            (self.mink, self.maxk) = self.experiment["ksweep"]
        self.bestk = 0
        self.delta = 0
        # Default initial size of sketch array for each node
        RANGEK = max(100, self.maxk+2)
        self.ksketches = [None] * RANGEK
        self.assign_progeny()
        self.fastas = [f.fastas[0] for f in self.progeny]
        self.ngen = len(self.progeny)
        

    def __repr__(self):
        return f"{self.__class__.__name__}['{self.node_title}', k: {self.bestk}, delta: {self.delta}, ngen: {self.ngen}, children: {repr(self.children)} ]"
        
    def __lt__(self, other):
        # lt = less than
        return self.ngen < other.ngen
    
    def assign_progeny(self):
        '''If a node doesn't have progeny, it is it's own progeny'''
        if not self.progeny:
            self.progeny=[self]
            self.fastas=[self.node_title]
            self.node_title=os.path.splitext(os.path.basename(self.node_title))[0]
            
    def find_delta_helper(self, kval: int, direction=1):
        '''determine if there is a local maximum delta relative to the current k'''

        if self.experiment['tool'] == 'dashing' and kval > 32:
            raise ValueError("Exploratory k value is too high for dashing. Either something is amiss with your data or you need to be using --exact mode")
        # what if maxk is larger than the default ksketches size?
        elif kval > len(self.ksketches):
            self.ksketches.extend([None]* (kval-len(self.ksketches)+2))
        # make sure that all necessary ingredient sketches are available for the current k in question
        self.node_ksweep(mink=min(kval-direction,kval,kval+direction), maxk=max(kval-direction,kval,kval+direction))
        self.update_node(kval)
        
        if direction < 0:
            self.mink = kval
        else:
            self.maxk = kval
        if self.delta == 0:
            if self.experiment["verbose"]:
                print("delta is 0 post update_node")
        old_d = self.delta
        new_d = self.ksketches[kval].delta_pos
        
        if old_d <= new_d:
            self.speciesinfo.kstart = kval
            self.bestk = kval
            self.delta = new_d
            self.find_delta_helper(kval=kval+direction, direction=direction)
        return
    
    def find_delta(self, kval: int):
        '''search in both directions of provided kvalue to detect a local maximum'''
        self.find_delta_helper(kval=kval, direction=1)
        self.find_delta_helper(kval=kval, direction=-1)
        self.card=self.ksketches[self.bestk].card
        return
    
    def ksweep_update_node(self, mink, maxk):
        # This should return an sfp that can be used to fill a parallel command
        sfp=SketchFilePath(filenames=self.fastas, kval=0, speciesinfo=self.speciesinfo, experiment=self.experiment)
        # If maxk is greater than the size of ksketches, expand it
        if maxk > len(self.ksketches):
            print("MAXK is larger than the size if ksketches, extending ksketches")
            self.ksketches.extend([None]* (maxk-len(self.ksketches)+2))
        # this paths list will be appended and passed during SketchObj creation
        presketches=[] 
        # this sketchlist will be used for a bath cardinality check
        sketchlist=[]
        # kmc sketches will need a parent tmp directory
        tmpdir=tempfile.mkdtemp()
        # use mink/maxk provided unless ksweep in experiment object is different
        krange=[int(mink),int(maxk)]
        
        # check the ks already in the node sketches, don't do them if they are already there. If none are empty return nothing
        empty_ks = [str(i) for i in range(int(mink),int(maxk)+1) if self.ksketches[i] is None]
        static_empty_ks = empty_ks.copy()
        if len(empty_ks) == 0:
            return []
        # create output directories for all the ks we are about to batch
        for i in range(int(krange[0]),int(krange[1])+1):
            os.makedirs(sfp.dir.replace("{}",str(i)), exist_ok=True)
            sketchlist.append(sfp.full.replace("{}",str(i)))
            self.experiment["baseset"].add(sfp.base.replace("{}",str(i)))
            
        if self.ngen > 1:
            for i in range(len(self.children)):
                sketchlist = sketchlist + self.children[i].ksweep_update_node(mink=mink, maxk=maxk)
                presketches= presketches + [self.children[i].ksketches[0].sfp.full]
        # store a "sketchobj" at k=0 that holds the parallel command 
        if self.experiment["tool"] == "dashing":
            self.ksketches[0] = DashSketchObj(kval = 0, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)
        elif self.experiment["tool"] == "kmc":
            self.ksketches[0] = KMCSketchObj(kval = 0, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)
            tmpdir=tempfile.mkdtemp()
            for k in empty_ks:
                os.mkdir(os.path.join(tmpdir,"k"+str(k)))
        # make a list of ks that don't have cardinalities/etc. stored in the tree and then figure out what the paths to those sketches would be
        for k in static_empty_ks:
            sketch_loc=self.ksketches[0].sfp.full.replace("{}",str(k))
            update_sketch = False
            # first check if those sketches exist
            if self.ksketches[0].sketch_check(path=sketch_loc):
                empty_ks.remove(k)
                if (sketch_loc in self.speciesinfo.cardkey.keys()):
                    if  self.speciesinfo.cardkey[sketch_loc] is not None and float(self.speciesinfo.cardkey[sketch_loc]) > 0 :
                        update_sketch = True
            elif self.experiment["lowmem"]:
                if (sketch_loc in self.speciesinfo.cardkey.keys()) and float(self.speciesinfo.cardkey[sketch_loc]) > 0 and self.ngen > 1:
                    empty_ks.remove(k)
                    update_sketch = True
            if update_sketch:
                pass

        if len(empty_ks) == 0:
            return []
        if self.ngen < 2:
            sub_cmd = self.ksketches[0]._leaf_command(tmpdir=tmpdir)
        else:
            sub_cmd = self.ksketches[0]._union_command()
        cmdlist = ["parallel -j 95% '",sub_cmd,"' :::",  " ".join(empty_ks)]
        cmd = " ".join(cmdlist)
        if self.experiment["debug"]:
            print(cmd)
        elif self.experiment["verbose"]:
            lines=cmd.splitlines()
            lnum=len(lines)
            i=0
            for i in range(lnum):
                if i<= 2 or i>=(lnum-2):
                    print(lines[i])
                elif i == (lnum-3):
                    print("\n")
                else:
                    print("..", end="")
            
        subprocess.call(cmd, shell=True, stdout=None)
        shutil.rmtree(tmpdir)
        
        return sketchlist

    def update_node(self, kval):
        '''Populate the sketch object for the given k at the self node as well as all children of the node'''
        ekrange=(kval, kval)
        if self.experiment["ksweep"] is not None:
            ekrange=self.experiment["ksweep"]
            if kval not in range(ekrange[0],ekrange[1]+1):
                print(f"k={kval} is outside of ksweep range ",ekrange)
                return
        if not self.ksketches[kval]:
            #create sketch file path holding information relating to the sketch for that k 
            sfp = SketchFilePath(filenames=self.fastas, kval=kval, speciesinfo=self.speciesinfo, experiment=self.experiment)
            # if this isn't a leaf node then collect the sketches for unioning
            presketches=[]
            if self.ngen > 1:
                presketches=[]
                # update each of the child nodes 
                for i in range(len(self.children)):
                    self.children[i].update_node(kval)
                    presketches= presketches + [self.children[i].ksketches[kval].sketch]
            # create a sketch object dependent on whether the tool is dashing or kmc
            if self.experiment["tool"] == "dashing":
                self.ksketches[kval] = DashSketchObj(kval = kval, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)
            elif self.experiment["tool"] == "kmc":
                self.ksketches[kval] = KMCSketchObj(kval = kval, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)
        return

    def node_ksweep(self, mink, maxk):
        '''Sketch all of the ks for the node (and its decendent nodes)between mink and maxk (even when they weren't needed to calculate delta'''
        sketchlist = self.ksweep_update_node(mink=mink, maxk=maxk)
        multicard = self.ksketches[0].card_command(sketchlist)
        self.ksketches[0].individual_card(cmd=multicard)

        # This loop is necessary to populate properties of the sketch objects
        for kval in range(mink, maxk+1):
            if self.ksketches[kval] is None:
                self.update_node(kval)

        self.mink = mink
        self.maxk = maxk
        return 

    def summarize(self, mink:int=0, maxk:int=0, ordering_number=0):
        '''create a dataframe of all "possible" delta values that were examined during creation of the node for use in plotting -- this isn't actually summarizing, so the function should be renamed'''
        nodevals=[]
        
        for kval in range(mink, maxk+1):
            linedict = {"ngen": self.ngen, "kval": kval, "card": self.ksketches[kval].card, "delta_pos": self.ksketches[kval].delta_pos, "title": self.node_title, "command": self.ksketches[kval].cmd, "ordering":ordering_number}
            nodevals.append(linedict)
        return nodevals 