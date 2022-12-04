#import parallel
#import random
from ast import Str
import sys
import os
import pickle
#import hashlib
#import warnings
import re
import subprocess
import csv
from multiprocessing import Process, Pool
from random import sample
from numpy import unique
import pandas as pd
from species_specifics import SpeciesSpecifics
from sketch_classes import SketchFilePath
from sketch_classes import SketchObj

#This is assuming that the command for dashing has been aliased
DASHINGLOC="dashing"
RANGEK=100

def retrieve_fasta_files(inputdir, full=True):
    '''return a list of all fasta files in a directory accounting for all the possible extensions'''
    reg_compile = re.compile(inputdir + "/*\.(fa.gz|fasta.gz|fna.gz|fasta|fa)")
    fastas = [fasta for fasta in os.listdir(inputdir) if reg_compile]
    if full:
        fastas=[os.path.join(inputdir,fasta) for fasta in fastas]
    return fastas

def random_orderings(length, norder, preexist=set()):
    '''create norder NEW unique random orderings of numbers 0 through length in addition to those in provided preexisting set of orderings'''
    rlist = list(range(length))
    prelen=len(preexist)
    newset=preexist.copy()
    totalneed = norder + prelen
    while len(newset) < totalneed:
        needed=totalneed - len(newset)
        #union the pre-existing set of orderings with remaining number of orderings required
        newset.update(set([tuple(sample(rlist,length)) for i in range(needed)]))
    return newset

def _update_helper(node, specinfo, kval):
    node.find_delta(specinfo, kval)
    return node.ksketches[kval].sketch

class DeltaTreeNode:
    ''' A node in a Delta tree. 
        node_title = name of input file or composite of inputfiles
        children = the nodes that are this nodes children
        progeny = list of leaf nodes decended from the node
        ksweep = should ks 1-100 be explored even after delta is found?
        
        '''
    def __init__(self, node_title: str, children: list, progeny=None ):
        self.node_title = node_title
        self.progeny = progeny
        # self.left = left
        # self.right = right
        self.children=children
        self.mink = 0
        self.maxk = 0
        self.bestk = 0
        self.delta = 0
        self.ksketches = [None] * RANGEK
        self.assign_progeny()
        self.fastas = [f.node_title for f in self.progeny]
        self.ngen = len(self.progeny)
        #print(self.progeny)
        #self.speciesinfo=speciesinfo

    def __repr__(self):
        return f"['{self.node_title}', k: {self.bestk}, delta: {self.delta}, ngen: {self.ngen} ]"
        
    def __lt__(self, other):
        # lt = less than
        return self.ngen < other.ngen
    # def is_leaf(self):
    #     if self.right or self.left:
    #         return False
    #     else:
    #         return True
    
    def assign_progeny(self):
        '''something'''
        if not self.progeny:
            self.progeny=[self]
    
            
#    def batch_update_card(self, card0, cards):
#        cmdlist = [DASHINGLOC,"card --presketched -p10"] +  card0
#        cmd = " ".join(cmdlist)
#        card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout.readlines()
#        for card in csv.DictReader(card_lines, delimiter='\t'):
#            cards[card['#Path']] = card['Size (est.)']
#        return cards
    def find_delta_helper(self, speciesinfo: SpeciesSpecifics, kval: int, direction=1):
        '''determine if there is a local maximum delta relative to the current k'''
        self.update_node(speciesinfo, kval)
        if direction < 0:
            self.mink = kval
        else:
            self.maxk = kval
        old_d = self.delta
        new_d = self.ksketches[kval].delta_pos
        if old_d < new_d:
            speciesinfo.kstart = kval
            self.bestk = kval
            self.delta = new_d
            self.find_delta_helper(speciesinfo=speciesinfo, kval=kval+direction, direction=direction)
        return
    
    def find_delta(self, speciesinfo: SpeciesSpecifics, kval: int):
        '''search in both directions of provided kvalue to detect a local maximum'''
        self.find_delta_helper(speciesinfo=speciesinfo, kval=kval, direction=1)
        #TODO maybe this one should reference the species kstart
        self.find_delta_helper(speciesinfo=speciesinfo, kval=kval, direction=-1)
        return

    def update_node(self, speciesinfo, kval, parallel=False):
        '''Populate the sketch object for the given k at the self node as well as all children of the node'''
        if self.ksketches[kval] is None:
            #create sketch file path holding information relating to the sketch for that k 
            sfp = SketchFilePath(filenames=self.fastas, kval=kval, speciesinfo=speciesinfo)
            # if this isn't a leaf node then collect the sketches for unioning
            if self.ngen > 1:
                presketches=[]
                # update each of the child nodes 
                # TODO can be done in parallel
                if parallel:
                    pool2 = Pool(processes=4)
                    results = [pool2.apply(_update_helper, args=(child, speciesinfo, kval, )) for child in self.children]
                    presketches=presketches + results
                else:
                    for i in range(len(self.children)):
                        self.children[i].update_node(speciesinfo, kval)
                        presketches= presketches + [self.children[i].ksketches[kval].sketch]

                #print(presketches)
                
                self.ksketches[kval] = SketchObj(kval = kval, sfp = sfp, speciesinfo=speciesinfo, presketches=presketches)

            # if this is a leaf node, sketch from fasta    
            elif self.ngen == 1:
                #print("Inside ngen=1 of update_node")
                self.ksketches[kval] = SketchObj(kval = kval, sfp = sfp,  speciesinfo=speciesinfo)
                self.update_card(speciesinfo)
            else:
                raise ValueError("For some reason you are trying to sketch ngen that is not >= 1. Something is amiss.")
                
        self.update_card(speciesinfo)
        return

    def ksweep(self, speciesinfo: SpeciesSpecifics, kmin: int, kmax: int):
        '''Sketch all of the ks for the node (and its decendent nodes)between kmin and kmax (even when they weren't needed to calculate delta'''
        if kmax > len(self.ksketches):
            self.ksketches = self.ksketches.extend([None]* (kmax-len(self.ksketches)))
        for kval in range(kmin, kmax+1):
            if not self.ksketches[kval]:
                self.update_node(speciesinfo, kval)
        

    def plot_df(self):
        '''create a dataframe of all "possible" delta values that were examined during creation of the node for use in plotting'''
        nodevals=[]
        for kval in range(len(self.ksketches)):
            if self.ksketches[kval]:
                linelist=[self.ngen,kval,self.ksketches[kval].delta_pos, self.node_title]
                nodevals.append(linelist)
        ndf=pd.DataFrame(nodevals,columns=["ngenomes","kval","delta_pos", "title"])
        return ndf

    def update_card(self, speciesinfo):
        '''retrieve/calculate cardinality for any sketches in the node that lack it'''
        sketches = [sketch for sketch in self.ksketches if sketch is not None]
        cardlist = [sketch.sketch for sketch in sketches if sketch.card == 0]
        if len(cardlist) > 0:
            cmdlist = [DASHINGLOC,"card --presketched -p10"] +  cardlist
            cmd = " ".join(cmdlist)
            print(cmd)
            card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout.readlines()
            for card in csv.DictReader(card_lines, delimiter='\t'):
                speciesinfo.cardkey[card['#Path']] = card['Size (est.)']
        #print(self._speciesinfo.cardkey)
        for sketch in self.ksketches:
            if sketch is not None:
                if sketch.sketch in sketches:
                    sketch.check_cardinality(speciesinfo)
                    sketch.delta_pos=sketch.card/sketch.kval

# TODO add argument for dashing v kmc --> estimate=True
class DeltaTree:
    ''' Delta tree data structure. '''
    def __init__(self, fasta_files, speciesinfo, nchildren=2):
        self.hashkey = speciesinfo.hashkey
        #self.files = self.fasta_files(speciesinfo.inputdir)
        #self.codebook = {}
        #self._code_lengths = []
        self._symbols = []
        #self._code_words = []
        self.mink=0
        self.maxk=0
        
        self.kstart=speciesinfo.kstart
        self.registers=speciesinfo.registers
        print(nchildren)
        self._build_tree(fasta_files, speciesinfo, nchildren)
        self.fill_tree(speciesinfo)
        self.ngen = len(fasta_files)
        self.root=self._dt[-1]
        self.delta = self.root_delta()
        self.fastas = fasta_files
        speciesinfo.kstart = self.root_k()
        self.speciesinfo=speciesinfo
        self.speciesinfo.save_hashkey()
        #self.card0 = speciesinfo.card0
        #self.cardkey=speciesinfo.cardkey
    def __sub__(self, other):
        # sub = subtraction
        print("Larger Tree Delta: ", self.delta)
        print("Subtree Delta: ", other.delta)
        print("subtraction result: ", self.delta - other.delta)
        return self.delta - other.delta

        
    def __repr__(self) -> None:
        ''' Traverse the DeltaTree in a depth-first way.'''
        root = self._dt[-1]
        print(root)
        def _print_tree_recursive(node) -> None:
            if node.children:
                nchild=len(node.children)
                print("NUMBER OF CHILDREN",nchild)
                for i in range(nchild):
                    
                    print(i)
                    n=node.children[i]
                    print("Child {}:".format(i),n)
                    _print_tree_recursive(n)
        _print_tree_recursive(root)

    ## NOTE batch_update_card not in use
    def batch_update_card(self, speciesinfo, debug=False):
        '''Update the cardinality dictionary for any sketches which have been added to the card0 list in speciesinfo '''
        if len(speciesinfo.card0) > 0:
            # TODO check for kmc v dashing
            cmdlist = [DASHINGLOC,"card --presketched -p10"] +  speciesinfo.card0
            cmd = " ".join(cmdlist)
            if debug:
                print(cmd)
            card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout.readlines()
            for card in csv.DictReader(card_lines, delimiter='\t'):
                speciesinfo.cardkey[card['#Path']] = card['Size (est.)']
            #print(speciesinfo.cardkey)
            for node in self._dt:
                for sketch in node.ksketches:
                    if sketch:
                        if sketch.sketch in speciesinfo.card0:
                            sketch.card=speciesinfo.check_cardinality(sketch.sketch)
                            sketch.delta_pos=sketch.card/sketch.kval
            speciesinfo.card0 = []
    
    def root_delta(self):
        '''retrieve the root node of the tree'''
        root=self._dt[-1]
        return root.delta
    def root_k(self):
        '''retrieve the argmax k for the root node'''
        root=self._dt[-1]
        return root.bestk


    def _build_tree(self, symbol: list, speciesinfo: SpeciesSpecifics, nchildren: int, parallel=False) -> None:
        '''
        Build a DeltaTree. The depth first nodes will have the provided number of children until there are only k<n input fastas left. A python list of nodes is returned with pointers to child nodes where applicable.

        When inserting a new node between node[3] and node[4]:

        * Python list:
            ```
            list = list[:3+1] + [new_node] + list[3+1:]
            ```

        Inputs:
            - symbol: a list of fasta files (str)
            - speciesinfo: SpeciesSpecifics object containing information specific to the overall species
            - nchildren: number of children of each node in the tree (or at least as many nodes as it works for)
        '''
        # create leaf nodes for all the provided fastas
        inputs = [
            DeltaTreeNode(
                node_title=s, children=None
            ) for s in symbol]
        inputs.sort()
        ## TODO: parallelize right here
        #for n in inputs:
        #    n.find_delta(speciesinfo, speciesinfo.kstart)
        #inputs = [n.find_delta(speciesinfo, self.registers, speciesinfo.kstart) for n in inputs]
        # 
        if parallel:
            pool = Pool(processes=4)
            results = [pool.apply(_update_helper, args=(dnode, speciesinfo,speciesinfo.kstart,  )) for dnode in inputs]
        else:
            for n in inputs:
                
                n.find_delta(speciesinfo, speciesinfo.kstart)
        
        
        self._dt = inputs
        idx_insert = 0
        idx_current = 0
        while idx_current != len(self._dt) - 1:
            increment=nchildren-1
            children=self._dt[idx_current:idx_current+nchildren]
            progeny=[p.progeny for p in children]
            #flatten the progeny list
            progeny=[item for sublist in progeny for item in sublist]
            child_titles=[c.node_title for c in children]
            new_node = DeltaTreeNode(
                node_title="_".join(child_titles),
                children = children,
                # children = [self._dt[idx_current],self._dt[idx_current+1]]
                progeny=progeny
            )
            new_node.find_delta(speciesinfo=speciesinfo, kval=speciesinfo.kstart)
            
            while idx_insert < len(self._dt)-increment and self._dt[idx_insert+increment].ngen <= new_node.ngen:
                
                idx_insert += increment
            
            self._dt = self._dt[:idx_insert+increment] + [new_node] + self._dt[idx_insert+increment:]
                
            idx_current += nchildren
            if idx_insert + increment > len(self._dt)-1:
                nchildren = len(self._dt) - idx_current
                #print("NCHILDREN:", nchildren)

        #self.batch_update_card()
        #self.compute_code()
    

    def print_list(self) -> None:
        nodes = []
        for i, node in enumerate(self._dt):
            nodes.append(f'\'{node.node_title}\'({node.ngen}\'({" ".join([i.node_title for i in node.progeny])})')
        print(' -> '.join(nodes))

    def fill_tree(self, speciesinfo):
        '''Starting at the root make sure that all nodes in the tree contain the sketches for the argmax ks for every node as well as 2 less than the minimum and 2 greater than the maximum'''
        root = self._dt[-1]
        bestks = list(unique([n.bestk for n in self._dt]))
        bestks = [k for k in bestks if k!=0  ]
        #print(bestks)
        bestks.sort()
        print("Best ks before padding:", bestks)
        bestks = bestks + [bestks[0]-1] + [bestks[0]-2] + [bestks[-1]+1] + [bestks[-1]+2]
        print("Best ks after padding:", bestks)
        for k in bestks:
            root.update_node( speciesinfo, k)
        
        
    def save(self, outloc: str, tag: str, label=None):
        '''Save the delta tree for future retrieval
        '''
        if label is None:
            label = ""
        else: label = "_"+label
        filepath=os.path.join(outloc,  tag + label + "_" + str(self.ngen) + '_dtree.pickle')
        with open(filepath,"wb") as f:
            pickle.dump(obj=self, file=f)
        print("Tree Pickle saved to: "+filepath)
        return filepath

  
    def delta_pos(self):
        ''' Traverse the DeltaTree to return a dataframe with all possible delta values.'''
        root = self._dt[-1]
        #print(root)
        def _delta_pos_recursive(node):
            tmplist=[node.plot_df()]
            
            if node.children:
                nchild=len(node.children)
                for i in range(nchild):
                    #print(i)
                    n=node.children[i]
                    tmplist.extend(_delta_pos_recursive(n))
            return tmplist
        
        dflist= _delta_pos_recursive(root)
        return pd.concat(dflist)

    def find_delta_delta(self, fasta_subset: list):
        '''Provided a list of fastas in a subset, find the delta-delta values between the whole spider and a spider without the provided fastas'''
        majord = self.delta
        majork = self.root_k()
        # create list of fastas that are in the original spider that are not in the subset provided
        fastas = [f for f in self.fastas if f not in fasta_subset]
        if len(fastas) == 0:
            fastas = [f for f in self.fastas if os.path.basename(f) not in fasta_subset]
        small_spider = DeltaSpider(fasta_files=fastas, speciesinfo=self.speciesinfo)
        print("Full Tree Delta: ", self.delta)
        print("Subtree Delta: ", small_spider.delta)
        return self - small_spider
   
    def ksweep(self, kmin=1, kmax=RANGEK):
        self.root.ksweep(self.speciesinfo, kmin=kmin, kmax=kmax)


class DeltaSpider(DeltaTree):
    '''Create a structure with all single sketches in terminal nodes tied to a single union node for all of them'''
    def __init__(self, fasta_files, speciesinfo):
        print("I made it to spider")
        nchildren=len(fasta_files)
        super().__init__(fasta_files, speciesinfo, nchildren=nchildren)


    def orderings_list(self, ordering_file=None, flist_loc=None, count=None):
        '''create or retrieve a series of random orderings of fasta sketches. return also the expected "sorted" array of the files. A subset of the fastas in the tree can be provided by name (in a file). The ordering of this file will be used when count=1 and the list is provided.'''
        fastas=self.fastas
        fastas.sort()

        # If a fasta file list is provided, subset the fastas from the species directory to only use the intersection
        if flist_loc:
            with open(flist_loc) as file:
                fsublist = [line.strip() for line in file]
            fastas = [f for f in fastas if f in fsublist]
            #order fastas as given in file
            fastas = [f for f in fsublist if f in fastas]
        # if count is one "sorted" ordering is returned with the reference list
        if count == 1:
            return [fastas,[[i for i in range(len(fastas))]]]
        
        # orderings are handled as sets to prevent duplication
        orderings=set()
        default_ordering=os.path.join(self.speciesinfo.sketchdir, self.speciesinfo.tag + "_"+ str(len(fastas))+"orderings.pickle")
        # if no ordering file is provided the default location is used
        if not ordering_file:
            ordering_file=default_ordering
        # if the ordering file exists then read the orderings
        if os.path.exists(ordering_file):
            with open(ordering_file,'rb') as f:
                orderings=pickle.load(f)
            # if count was not provided then just return the orderings that are already there
            if not count:
                return [fastas, list(orderings)]
            # if the count is lte to the number of orderings in the file, take the first count number of orderings
            if count <= len(orderings):
                return [fastas, list(orderings)[:count]]

        # if there is no ordering file at the location, time to make some    
        else:
            # if count is not provided, how will we know how many to make??
            if not count:
                raise ValueError("You must provide a value for count when there is no default ordering file")
        
        orderings = random_orderings(length=len(fastas), norder=count-len(orderings), preexist=orderings)
        # save the orderings for use next run of species 
        with open(ordering_file,"wb") as f:
            pickle.dump(orderings, f)
        return fastas, list(orderings)

    def progressive_wrapper(self, flist_loc=None, count=30, ordering_file=None,step=1):
        fastas, orderings = self.orderings_list( ordering_file=ordering_file, flist_loc=flist_loc, count=count)

        return self.progressive_union(flist=fastas, orderings=orderings, step=step)

    def progressive_union(self, flist, orderings, step):
        '''create (or use if provided) a series of random orderings to use when adding the individual fasta sketches to a union. Outputs a table with the delta values and associated ks at each stage'''

        #speciesinfo = SpeciesSpecifics(tag=tag, genomedir=genomedir, sketchdir=sketchdir, kstart=kstart)
        
        #fastas = retrieve_fasta_files(speciesinfo.inputdir)
        #fastas=self.fastas
        # fastas, orderings = self.orderings_list( ordering_file=ordering_file, flist_loc=flist_loc, count=count)
        # speciesinfo=self.speciesinfo
        

        # create a sketch of the full union of the fastas
        print(flist)
        smain = DeltaSpider(fasta_files=flist, speciesinfo=self.speciesinfo)
        results=[]
        for i in range(0,len(orderings)):
            oresults=smain.sketch_ordering(orderings[i], step=step)
            odf=pd.DataFrame(oresults,columns=["ngenomes","kval","delta"])
            odf['ordering']=i+1
            results.append(odf)
            self.speciesinfo.save_hashkey()
        
        #rdf=pd.DataFrame(results,columns=["k","delta"])
        return pd.concat(results)


    def sketch_ordering(self, ordering, step=1):
        '''Provided an ordering for the fastas in a tree, create sketches of the subsets within that ordering and report the deltas in a dataframe'''
        flen=len(ordering)
        print(flen)
        output=[]
        for i in range(1,flen+1):
            if i % step == 0:
                sublist=[self.fastas[j] for j in ordering[:i]]
                ospider=DeltaSpider(fasta_files=sublist, speciesinfo=self.speciesinfo)
                output.append([i, ospider.root_k(), ospider.delta])
        return output


# class ProgressiveUnion:
#     '''something'''
#     def __init__(self, fasta_files, speciesinfo, nchildren=2):
#         self.hashkey = speciesinfo.hashkey
# TODO: add arg for dashing vs kmc
def create_delta_tree(tag: str, genomedir: str, sketchdir: str, kstart: int, nchildren=None, registers=None, flist_loc=None, canonicalize=True):
    '''Given a species tag and a starting k value retrieve a list of fasta files to create a tree with the single fasta sketches populating the leaf nodes and the higher level nodes populated by unions
    tag = species tag
    genomedir = parent directory of species subdirectory
    sketchdir = parent directory where species directory for output sketches should be created
    kstart = starting k to use while searching for delta
    nchildren = number of children that nodes should have (until they can't)
    registers = number of registers to use when sketching
    flist_loc = file containing list of subset of fasta files to use from species directory (IN FUTURE maybe list of fastas with loc?)'''
    # create a SpeciesSpecifics object that will tell us where the input files can be found and keep track of where the output files should be written
    
    speciesinfo = SpeciesSpecifics(tag=tag, genomedir=genomedir, sketchdir=sketchdir, kstart=kstart, registers=registers, flist_loc=flist_loc, canonicalize=canonicalize)
    print("HERE")
    #inputdir = speciesinfo.inputdir
    if flist_loc:
        with open(flist_loc) as file:
            fastas = [line.strip() for line in file]
    # right now we expect that if we are provided with a genome directory AND a file list the file list will only contain the basenames
    elif os.path.exists(speciesinfo.inputdir):
        print("I think the input directory exists")
        fastas = retrieve_fasta_files(speciesinfo.inputdir, full=True)
        #fastas = [f for f in allfastas if f in fastas]
    # If a fasta file list is provided subset the fastas from the species directory to only use the intersection
    # if flist_loc:
    #     with open(flist_loc) as file:
    #         fsublist = [line.strip() for line in file]
    #     fastas = [f for f in fastas if f in fsublist]
    fastas.sort()
    if nchildren:
        dtree = DeltaTree(fasta_files=fastas,speciesinfo=speciesinfo, nchildren=nchildren)
    else:
        dtree = DeltaSpider(fasta_files=fastas,speciesinfo=speciesinfo)
    # Save the cardinality keys as well as the hashkey for the next run of the species
    speciesinfo.save_cardkey()
    speciesinfo.save_hashkey()
    print(dtree)
    return dtree