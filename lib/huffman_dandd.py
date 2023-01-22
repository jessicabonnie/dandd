#import parallel
#import random
#from ast import Str
# import sys
import os
import pickle
#import hashlib
#import warnings
#import re
import subprocess
import csv
# from multiprocessing import Process, Pool
from random import sample, shuffle
# from numpy import unique
# import pandas as pd
from species_specifics import SpeciesSpecifics
from sketch_classes import SketchFilePath
from sketch_classes import SketchObj, DashSketchObj, KMCSketchObj
# import tabulate
from typing import List, Dict, Set, Tuple, NamedTuple
from string import ascii_uppercase
from dandd_cmd import write_listdict_to_csv
from math import factorial
from itertools import permutations


#This is assuming that the command for dashing has been aliased
DASHINGLOC="dashing" #"/scratch16/blangme2/jessica/lib/dashing/dashing"
RANGEK=100

# def permutations(elements):
#     '''Generate all permutations of a list of elements or a string'''
#     if len(elements) <= 1:
#         yield elements
#         return
#     for perm in permutations(elements[1:]):
#         for i in range(len(elements)):
#             # nb elements[0:1] works in both string and list contexts
#             yield perm[:i] + elements[0:1] + perm[i:]

def permute(length, norder, preexist=set())-> Set[Tuple[int]]:
    '''Create a set of ordering tuples'''
    fact=factorial(length)
    newset=preexist.copy()
    norder = min(norder, fact)
    print(f"{norder} permutations will be produced.")
    # if the number of permutations is low enough, generate all of them
    if fact < 4000:
        print("fact less than 4000")
        # create a randomized list of all permutations
        newpermute=list(permutations(range(length)))
        shuffle(newpermute)
    else:
        newpermute=list()
        while len(newpermute) < norder:
            # tmpset=set([tuple(sample(list(range(length)),length)) for _ in range(norder)])
            newpermute.extend([tuple(sample(list(range(length)),length)) for _ in range(norder)])
            newpermute=list(set(newpermute))
    index = 0
    while len(newset) < norder:
        newset.update([newpermute[index]])
        index+=1
    return newset

# def random_orderings(length, norder, preexist=set())->Set[Tuple[int]]:
#     '''create norder NEW unique random orderings of numbers 0 through length in addition to those in provided preexisting set of orderings'''
#     rlist = list(range(length))
#     prelen=len(preexist)
#     newset=preexist.copy()
#     totalneed = norder + prelen
#     maxpermute=factorial(length)
#     if totalneed > maxpermute:
#         ##TODO: check to see if norder needs to be changed? Maybe do this sooner?
#         totalneed=maxpermute
#         print(f"You asked for more permutations than are possible with the set size. Number of permutations is being reset to: {maxpermute}")
#     while len(newset) < totalneed:
#         needed=totalneed - len(newset)
#         #union the pre-existing set of orderings with remaining number of orderings required
#         newset.update(set([tuple(sample(rlist,length)) for i in range(needed)]))
#     return newset

# def _update_helper(node, kval):
#     node.find_delta(kval=kval)
#     return node.ksketches[kval].sketch

class DeltaTreeNode:
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
        self.bestk = 0
        self.delta = 0
        self.ksketches = [None] * RANGEK
        self.assign_progeny()
        self.fastas = [f.node_title for f in self.progeny]
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
    

    def find_delta_helper(self, kval: int, direction=1):
        '''determine if there is a local maximum delta relative to the current k'''
        # make sure that all necessary ingredient sketches are available for the current k in question
        self.update_node(kval)
        if direction < 0:
            self.mink = kval
        else:
            self.maxk = kval
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
        # seek it after a small jump also
        #self.find_delta_helper(kval = self.bestk + 2, direction=1)
        #self.find_delta_helper(kval = self.bestk + 4, direction=1)
        #TODO maybe this one should reference the species kstart
        self.find_delta_helper(kval=kval, direction=-1)
        self.card=self.ksketches[self.bestk].card
        return

    def update_node(self, kval, parallel=False):
        '''Populate the sketch object for the given k at the self node as well as all children of the node'''
        if not self.ksketches[kval]:
            #create sketch file path holding information relating to the sketch for that k 
            sfp = SketchFilePath(filenames=self.fastas, kval=kval, speciesinfo=self.speciesinfo, experiment=self.experiment)
            # if this isn't a leaf node then collect the sketches for unioning
            presketches=None
            if self.ngen > 1:
                presketches=[]
                # update each of the child nodes 
                # TODO can be done in parallel
                if parallel:
                    pass
                    # pool2 = Pool(processes=4)
                    # results = [pool2.apply(_update_helper, args=(child, self.speciesinfo, kval, )) for child in self.children]
                    # presketches=presketches + results
                else:
                    for i in range(len(self.children)):
                        self.children[i].update_node(kval)
                        presketches= presketches + [self.children[i].ksketches[kval].sketch]

            if self.experiment["tool"] == "dashing":
                self.ksketches[kval] = DashSketchObj(kval = kval, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)
            elif self.experiment["tool"] == "kmc":
                # TODO: create tmp directory for this obj to use here?
                self.ksketches[kval] = KMCSketchObj(kval = kval, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)

            # # if this is a leaf node, sketch from fasta    
            # elif self.ngen == 1:
            #     #print("Inside ngen=1 of update_node")
            #     self.ksketches[kval] = SketchObj(kval = kval, sfp = sfp,  speciesinfo=self.speciesinfo, experiment=self.experiment)
            #     self.update_card()
            # else:
            #     raise ValueError("For some reason you are trying to sketch ngen that is not >= 1. Something is amiss.")
                
        #self.update_card()
        return

    def node_ksweep(self, mink: int, maxk: int):
        '''Sketch all of the ks for the node (and its decendent nodes)between mink and maxk (even when they weren't needed to calculate delta'''
        print(len(self.ksketches))
        maxk=int(maxk)
        mink=int(mink)
        if maxk > len(self.ksketches):
            self.ksketches = self.ksketches + [None]* (maxk-len(self.ksketches))
            #self.ksketches.extend([None]* (maxk-len(self.ksketches)))
        for kval in range(mink, maxk+1):
            if not self.ksketches[kval]:
                self.update_node(kval) 
        self.mink=mink
        self.maxk=maxk    

    def summarize(self, mink:int=0, maxk:int=0):
        '''create a dataframe of all "possible" delta values that were examined during creation of the node for use in plotting -- this isn't actually summarizing, so the function should be renamed'''
        nodevals=[]
        if mink == 0:
            mink=self.mink
        if maxk == 0:
            maxk=self.maxk
        labels=[ascii_uppercase[i]+ ascii_uppercase[j] for i in range(26) for j in range(26)]

        for kval in range(mink, maxk+1):
            if self.ksketches[kval]:
                linedict = {"ngen": self.ngen, "kval": kval, "card": self.ksketches[kval].card, "delta_pos": self.ksketches[kval].delta_pos, "title": self.node_title}
                for index, value in enumerate(self.fastas):
                    linedict[f"step{index}"] = value
                    # linedict[labels[index]] = value
                
                # if len(self.fastas) <= 2:
                #     linedict.update({"A": self.fastas[0]})
                # if len(self.fastas) == 2:
                #     linedict.update({"B": self.fastas[1]})
                
                # linelist=[self.ngen,kval, self.ksketches[kval].card, self.ksketches[kval].delta_pos, self.node_title]
                nodevals.append(linedict)
        # ndf=pd.DataFrame(nodevals,columns=["ngenomes","kval","card","delta_pos", "title"])
        # return ndf
        return nodevals

    def update_card(self):
        '''retrieve/calculate cardinality for any sketches in the node that lack it. NOTE: this is not in use. also it is more sketch related. TODO: refactor to use commands from SketchObj classes rather than using dashing '''
        sketches = [sketch for sketch in self.ksketches if sketch is not None]
        cardlist = [sketch.sketch for sketch in sketches if sketch.card == 0]
        if len(cardlist) > 0:
            cmdlist = [DASHINGLOC,"card --presketched -p10"] +  cardlist
            cmd = " ".join(cmdlist)
            print(cmd)
            card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout.readlines()
            for card in csv.DictReader(card_lines, delimiter='\t'):
                self.speciesinfo.cardkey[card['#Path']] = card['Size (est.)']
        for sketch in self.ksketches:
            if sketch is not None:
                if sketch.sketch in sketches:
                    sketch.check_cardinality()
                    sketch.delta_pos=sketch.card/sketch.kval

class DeltaTree:
    ''' Delta tree data structure. '''
    def __init__(self, fasta_files, speciesinfo, nchildren=2, experiment={'tool':'dashing', 'registers':20, 'canonicalize':True, 'debug':False, 'nthreads':10, 'baseset': set(), 'safety': False}, padding=True):
        # self.fastahex = speciesinfo.fastahex
        self.experiment=experiment
        self._symbols = []
        #self._code_words = []
        self.mink=0
        self.maxk=0
        
        self.kstart=speciesinfo.kstart
        self.speciesinfo=speciesinfo
        self._build_tree(fasta_files, nchildren)
        self.fill_tree(padding=padding)
        self.ngen = len(fasta_files)
        self.root=self._dt[-1]
        self.delta = self.root_delta()
        self.fastas = fasta_files
        speciesinfo.kstart = self.root_k()
        self.speciesinfo.save_references()
    def __sub__(self, other):
        # sub = subtraction
        print("Larger Tree Delta: ", self.delta)
        print("Subtree Delta: ", other.delta)
        print("subtraction result: ", self.delta - other.delta)
        return self.delta - other.delta
    def __repr__(self):
       """Return a string which when eval'ed will rebuild tree"""
       return '{}(FASTAS: {}, NODES: {})'.format(
                self.__class__.__name__,
                self.fastas,
                repr(self._dt[-1]))

    def print_tree(self):
        ''' Traverse the DeltaTree in a depth-first way.'''
        root = self._dt[-1]
        print(root)
        def _print_tree_recursive(node):
            if not node.children.is_empty():
                nchild=len(node.children)
                print("NUMBER OF CHILDREN",nchild)
                for i in range(nchild):
                    
                    print(i)
                    n=node.children[i]
                    print("Child {}:".format(i),n)
                    _print_tree_recursive(n)
        _print_tree_recursive(root)

    ## NOTE batch_update_card not in use
    def batch_update_card(self):
        '''Update the cardinality dictionary for any sketches which have been added to the card0 list in speciesinfo '''
        if len(self.speciesinfo.card0) > 0:
            # TODO check for kmc v dashing
            cmdlist = [DASHINGLOC,"card --presketched -p10"] +  self.speciesinfo.card0
            cmd = " ".join(cmdlist)
            if self.experiment['debug']:
                print(cmd)
            card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout.readlines()
            for card in csv.DictReader(card_lines, delimiter='\t'):
                self.speciesinfo.cardkey[card['#Path']] = card['Size (est.)']
            for node in self._dt:
                for sketch in node.ksketches:
                    if sketch:
                        if sketch.sketch in self.speciesinfo.card0:
                            sketch.card=self.speciesinfo.check_cardinality(sketch.sketch)
                            sketch.delta_pos=sketch.card/sketch.kval
            self.speciesinfo.card0 = []
    
    def root_delta(self):
        '''retrieve the root node of the tree'''
        root=self._dt[-1]
        return root.delta
    def root_k(self):
        '''retrieve the argmax k for the root node'''
        root=self._dt[-1]
        return root.bestk


    def _build_tree(self, symbol: list, nchildren: int, parallel=False) -> None:
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
                node_title=s, children=[], speciesinfo=self.speciesinfo, experiment=self.experiment, progeny=[]
            ) for s in symbol]
        inputs.sort()
        ## TODO: parallelize right here
        #for n in inputs:
        #    n.find_delta(speciesinfo, speciesinfo.kstart)
        #inputs = [n.find_delta(speciesinfo, self.registers, speciesinfo.kstart) for n in inputs]
        # 
        if parallel:
            # self.parallel_progeny_prep(kval)
            pass
            # pool = Pool(processes=4)
            # results = [pool.apply(_update_helper, args=(dnode, self.speciesinfo,self.speciesinfo.kstart,  )) for dnode in inputs]
        else:
            for n in inputs:
                
                n.find_delta(self.speciesinfo.kstart)
        
        
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
                node_title="_".join(child_titles), speciesinfo=self.speciesinfo,
                children = children,
                progeny=progeny,
                experiment=self.experiment
            )
            new_node.find_delta(kval=self.speciesinfo.kstart)
            
            while idx_insert < len(self._dt)-increment and self._dt[idx_insert+increment].ngen <= new_node.ngen:
                idx_insert += increment
            
            self._dt = self._dt[:idx_insert+increment] + [new_node] + self._dt[idx_insert+increment:]
                
            idx_current += nchildren
            if idx_insert + increment > len(self._dt)-1:
                nchildren = len(self._dt) - idx_current

        #self.batch_update_card()
        #self.compute_code()
    
    # def parallel_progeny_prep(self, kval:int):
    #     cmd = parallel_progeny_command(self.speciesinfo.sketchdir, kval, self.experiment)
    #     proc = subprocess.Popen(cmd,shell=True)
    #     output, errs=proc.communicate(input="\n".join(self.fastas).encode())

    def print_list(self) -> None:
        nodes = []
        for i, node in enumerate(self._dt):
            nodes.append(f'\'{node.node_title}\'({node.ngen}\'({" ".join([i.node_title for i in node.progeny])})')
        print(' -> '.join(nodes))

    def fill_tree(self, padding=False):
        '''Starting at the root make sure that all nodes in the tree contain the sketches for the argmax ks for every node as well as 2 less than the minimum and 2 greater than the maximum (IF padding argument is True)'''
        root = self._dt[-1]
        bestks = list(set([n.bestk for n in self._dt]))
        bestks = [k for k in bestks if k!=0  ]
        #print(bestks)
        bestks.sort()
        # print("Best ks before padding:", bestks)
        if padding:
            bestks = bestks + [bestks[0]-1] + [bestks[0]-2] + [bestks[-1]+1] + [bestks[-1]+2] + [bestks[-1]+3]
        # print("Best ks after padding:", bestks)
        for k in bestks:
            root.update_node(k)
    
    def leaf_nodes(self) -> List[DeltaTreeNode]:
        return [child for child in self._dt if child.ngen==1]

    def to_spider(self):
        '''Transform tree into a spider if it isn't'''
        #TODO: make this use the leaf nodes of self
        # dspider = DeltaSpider(fasta_files=self.fastas,speciesinfo=self.speciesinfo,experiment=self.experiment)
        # self._dt = self.leaf_nodes()

        children=self.leaf_nodes()
        progeny=[p.progeny for p in children]
        #flatten the progeny list
        progeny=[item for sublist in progeny for item in sublist]
        child_titles=[c.node_title for c in children]
        body_node = DeltaTreeNode(
            node_title="_".join(child_titles), speciesinfo=self.speciesinfo,
            children = children,
            progeny=progeny,
            experiment=self.experiment
            )
        body_node.find_delta(kval=self.speciesinfo.kstart)
        self._dt = children + [body_node]
        
    # def get_experiment_map(self, sketchinfo):
    #     return sketchinfo.fromkeys(self.experiment["baseset"])
    
    # def name_exmap(self, outdir, tag, label):
        # bases=list(self.experiment["baseset"])
        # print(bases)
        # explist=[self.speciesinfo.sketchinfo[item] for item in bases]
        # keys=explist[0].keys()
        # filepath=os.path.join(outdir,  tag + label + "_" + str(self.ngen) + "_" + self.experiment["tool"] + '_sketchdb.txt')
        # return filepath
        # with open(filepath, "w") as writer:
        #     dict_writer = csv.DictWriter(writer, fieldnames=keys)
        #     dict_writer.writeheader()
        #     dict_writer.writerows(explist)


    def make_prefix(self, tag: str, label="", outdir:str=None):
        if not outdir:
            outdir=os.path.curdir
        if not label == "":
            label = "_"+label
        fileprefix=os.path.join(outdir,  tag + label + "_" + str(self.ngen) + "_" + self.experiment["tool"] )
        return fileprefix

    def save(self, fileprefix:str):
        '''Save the delta tree for future retrieval
        '''
        # if not label == "":
        #     label = "_"+label
        # fileprefix=self.make_prefix(tag=tag, label=label, outdir=outdir)
        
        filepath=fileprefix + '_dtree.pickle'
        with open(filepath,"wb") as f:
            pickle.dump(obj=self, file=f)
        print("Tree Pickle saved to: "+filepath)
        #subset the fastahex map to output a human readable version containing info for the sketches relevant to the tree
        expmaploc=fileprefix + '_sketchdb.txt'
        explist=[self.speciesinfo.sketchinfo[item] for item in list(self.experiment["baseset"])]
        write_listdict_to_csv(outfile=expmaploc, listdict=explist)
        print(f"Output Sketch/DB mapping saved to {expmaploc}.")
        return filepath

  
    def summarize(self, mink=0, maxk=0):
        ''' Traverse the DeltaTree to return a dataframe with all possible delta values. -- this isn't actually summarizing, so the function should be renamed'''
        root = self._dt[-1]
        if mink == 0:
            mink=self.mink
            if self.mink > 4:
                mink=self.mink - 2
        if maxk == 0:
            maxk = self.maxk
            if self.maxk <= 30:
                maxk=self.maxk + 2
        self.ksweep(mink=mink, maxk=maxk)
        #print(root)
        def _delta_pos_recursive(node):
            tmplist=node.summarize()
            
            if node.children:
                nchild=len(node.children)
                for i in range(nchild):
                    #print(i)
                    n=node.children[i]
                    tmplist.extend(_delta_pos_recursive(n))
            return tmplist
        
        dictlist= _delta_pos_recursive(root)
        #print(dictlist)
        return dictlist

    def nodes_from_fastas(self, fasta_list):
        return [node for node in self.leaf_nodes() if node.fastas[0] in fasta_list]

    def find_delta_delta(self, fasta_subset: List[str]) -> float:
        '''Provided a list of fastas in a subset, find the delta-delta values between the whole spider and a spider without the provided fastas'''
        # majord = self.delta
        # majork = self.root_k()
        # create list of fastas that are in the original spider that are not in the subset provided --> i.e. the complement
        fastas = [f for f in self.fastas if f not in fasta_subset]

        # if len(fastas) == 0:
        #     fastas = [f for f in self.fastas if os.path.basename(f) not in fasta_subset]
        small_spider = SubSpider(leafnodes=self.nodes_from_fastas(fastas), speciesinfo=self.speciesinfo, experiment=self.experiment)
        print("Full Tree Delta: ", self.delta)
        print("Subtree Delta: ", small_spider.delta)
        return self - small_spider
   
    def ksweep(self, mink:int=0, maxk:int=0):
        if mink == 0:
            mink=int(self.mink)
        if maxk == 0:
            maxk = int(self.maxk)
        self.root.node_ksweep( mink=mink, maxk=maxk)

    def orderings_list(self, ordering_file=None, flist_loc=None, count=0)-> Tuple[List[str], List[Tuple[int]]]:
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
            return fastas,[tuple(i for i in range(len(fastas)))]
        
        # orderings are handled as sets to prevent duplication
        orderings=set()
        default_ordering=os.path.join(self.speciesinfo.sketchdir, self.speciesinfo.tag + "_"+ str(len(fastas))+"_orderings.pickle")
        # if no ordering file is provided the default location is used
        if not ordering_file:
            ordering_file=default_ordering
        # if the ordering file exists then read the orderings
        if os.path.exists(ordering_file):
            with open(ordering_file,'rb') as f:
                orderings=pickle.load(f)
            # if count was not provided then just return the orderings that are already there
            if count==0:
                return fastas, list(orderings)
            # if the count is lte to the number of orderings in the file, take the first count number of orderings
            if count <= len(orderings):
                return fastas, list(orderings)[:count]

        # if there is no ordering file at the location, time to make some    
        else:
            # if count is not provided, how will we know how many to make??
            if count<1:
                raise ValueError("You must provide a value for count when there is no default ordering file")
        
        orderings = permute(length=len(fastas), norder=count, preexist=orderings)
        # orderings = random_orderings(length=len(fastas), norder=count-len(orderings), preexist=orderings)
        # save the orderings for use next run of species 
        with open(ordering_file,"wb") as f:
            pickle.dump(orderings, f)
        return fastas, list(orderings)

    def progressive_wrapper(self, flist_loc=None, count=30, ordering_file=None,step=1, debug=False)-> List[dict]:
        fastas, orderings = self.orderings_list( ordering_file=ordering_file, flist_loc=flist_loc, count=count)

        return self.progressive_union(flist=fastas, orderings=orderings, step=step)

    def progressive_union(self, flist, orderings, step) -> List[dict]:
        '''create (or use if provided) a series of random orderings to use when adding the individual fasta sketches to a union. Outputs a table with the delta values and associated ks at each stage'''

        # TODO: Update experiment object if needed
        # create a sketch of the full union of the fastas
        smain = DeltaSpider(fasta_files=flist, speciesinfo=self.speciesinfo, experiment=self.experiment)
        results=[]
        summary=[]
        for i in range(0,len(orderings)):
            oresults, osummary =smain.sketch_ordering(orderings[i], number=i+1, step=step)
            for o in osummary:
                o["ordering"] = i+1
            results.extend(oresults)
            summary.extend(osummary)
            self.speciesinfo.save_references()
            self.speciesinfo.save_cardkey(tool=self.experiment["tool"])
        return results, summary

    def sketch_ordering(self, ordering, number, step=1)->List[dict]:
        '''Provided an ordering for the fastas in a tree, create sketches of the subsets within that ordering and report the deltas'''
        flen=len(ordering)
        output=[]
        summary=[]
        for i in range(1,flen+1):
            if i % step == 0:
                sublist=[self.fastas[j] for j in ordering[:i]]
                ospider=SubSpider(leafnodes=self.nodes_from_fastas(sublist), speciesinfo=self.speciesinfo, experiment=self.experiment)
                output.append({"ngen":i, "kval":ospider.root_k(), "delta": ospider.delta, "ordering": number, "fastas": sublist})
                summary.extend(ospider.summarize())

        return output, summary

    
    def pairwise_spiders(self, sublist=[], mink=0, maxk=0):
        '''Create values for K-Independent-Jaccard. '''
        # super_spider=self.to_spider()

        if len(sublist)==0:
            sublist=self.leaf_nodes()
        pairings=[[a, b] for idx, a in enumerate(sublist) for b in sublist[idx + 1:]]
        kij_results=[]
        j_results=[]
        for pair in pairings:

        # def pairwise_helper(pair):
            pspider=SubSpider(leafnodes=pair,speciesinfo=self.speciesinfo,experiment=self.experiment)
            pspider.ksweep(mink=mink, maxk=maxk)
            kij_results.append(pspider.kij_summarize())
            j_results.extend(pspider.jaccard_summarize(mink=mink, maxk=maxk))
            # pspider.ksweep(mink=mink, maxk=maxk)
            # childA=pspider._dt[0]
            # childB=pspider._dt[1]
            # outdict={"A":childA.fastas[0], "B":childB.fastas[0],
            # "Adelta":childA.delta, "Bdelta":childB.delta,
            # "Ak":childA.bestk, "Bk":childB.bestk,
            # "ABdelta":pspider.delta, "ABk":pspider.root.bestk}
            # outdict["KIJ"]=(outdict["Adelta"] + outdict["Bdelta"]-outdict["ABdelta"])/outdict["ABdelta"]
            
        # results=[pairwise_helper(pair) for pair in pairings]
        self.speciesinfo.save_references()
        self.speciesinfo.save_cardkey(tool=self.experiment["tool"])
        return kij_results, j_results
    

class SubSpider(DeltaTree):
    def __init__(self,leafnodes,speciesinfo,experiment):
        self.fastahex = speciesinfo.fastahex
        self.experiment=experiment
        self._symbols = []
        # self.mink=0
        # self.maxk=0
        self.kstart=speciesinfo.kstart
        self.speciesinfo=speciesinfo
        self._build_tree(leafnodes)
        self.root=self._dt[-1]
        self.fastas=self.root.fastas
        self.ngen = len(self.fastas)
        self.delta = self.root_delta()
        self.speciesinfo.save_references()
        #super().__new__(self)

    def _build_tree(self, leafnodes):
        children=leafnodes
        progeny=[p.progeny for p in children]
        #flatten the progeny list
        progeny=[item for sublist in progeny for item in sublist]
        child_titles=[c.node_title for c in children]
        body_node = DeltaTreeNode(
            node_title="_".join(child_titles), speciesinfo=self.speciesinfo,
            children = children,
            progeny=progeny,
            experiment=self.experiment
            )
        body_node.find_delta(kval=self.speciesinfo.kstart)
        self.mink=body_node.mink
        self.maxk=body_node.maxk
        self._dt = children + [body_node]

    def kij_summarize(self):
        '''Calculate k independent jaccard for the subspider'''
        ##TODO: Write this to handle more than 2 children?
        if len(self.fastas) != 2:
            raise ValueError("KIJ can only be calculated on spider/trees with 2 children")
        childA=self._dt[0]
        childB=self._dt[1]
        outdict={"A":childA.fastas[0], "B":childB.fastas[0],
            "Adelta":childA.delta, "Bdelta":childB.delta,
            "Ak":childA.bestk, "Bk":childB.bestk,
            "ABdelta":self.delta, "ABk":self.root.bestk}
        outdict["KIJ"]=(outdict["Adelta"] + outdict["Bdelta"]-outdict["ABdelta"])/outdict["ABdelta"]
        return outdict
    
    def jaccard_summarize(self, mink=0, maxk=0):
        '''Calculate jaccard distance for the subspider'''
        ##TODO: Write this to handle more than 2 children?
        if len(self.fastas) != 2:
            raise ValueError("KIJ can only be calculated on spider/trees with 2 children")
        childA=self._dt[0]
        childB=self._dt[1]
        jevals=[]
        outdict={"A":childA.fastas[0], "B":childB.fastas[0]}
        # self.ksweep(mink=mink,maxk=maxk)
        for k in range(mink,maxk+1):
            outdict.update({"kval":k , "Acard": childA.ksketches[k].card, "Bcard": childB.ksketches[k].card, "ABcard":self.root.card})
            outdict["jaccard"] = (outdict["Acard"] + outdict["Bcard"] - outdict["ABcard"])/outdict["ABcard"]
            jevals.append(outdict)
        return jevals

class DeltaSpider(DeltaTree):
    '''Create a structure with all single sketches in terminal nodes tied to a single union node for all of them'''
    def __init__(self, fasta_files, speciesinfo, experiment, padding=False):
        # print("I made it to spider")
        nchildren=len(fasta_files)
        super().__init__(fasta_files=fasta_files, speciesinfo=speciesinfo, experiment=experiment, nchildren=nchildren, padding=padding)
    def __init2__(self, tree:DeltaTree):
        raise NotImplementedError("initialization of spider by tree not yet implemented")
    #TODO: add a function to the superclass that adds union nodes to the _dt and then create one here that adds just the spider body. That way the function can be passed the list of childnodes and then repurposed to instantiate a spider using only child nodes
    # def _build_tree(self.leaf_nodes=[DeltaTreeNode]):

    ## TODO: init function receives a delta tree and creates a spider out of it's nodes without creating new child nodes

# class ProgressiveUnion:
#     '''Unimplemented Progressive Union Class to hold functions and objects relating to that capability'''
#     def __init__(self, deltatree: DeltaTree, orderings: list):
#         self.orderings=[]
#         self.dtree=deltatree
#         self.orderings=orderings

#         raise NotImplementedError


def create_delta_tree(tag: str, genomedir: str, sketchdir: str, kstart: int, nchildren=None, registers=0, flist_loc=None, canonicalize=True, tool='dashing', debug=False, nthreads=10, safety=False):
    '''Given a species tag and a starting k value retrieve a list of fasta files to create a tree with the single fasta sketches populating the leaf nodes and the higher level nodes populated by unions
    tag = species tag
    genomedir = parent directory of species subdirectory
    sketchdir = parent directory where output sketches should be created
    kstart = starting k to use while searching for delta
    nchildren = number of children that nodes should have (until they can't)
    registers = number of registers to use when sketching
    flist_loc = file containing list of subset of fasta files to use from species directory (IN FUTURE maybe list of fastas with loc?)
    canonicalize = T/F indicating whether kmers should be canonicalized
    tool = string indicating which tool to use for kmer cardinality
    choices=["dashing","kmc"] '''
    # create an experiment dictionary for values that are needed at multiple levels that are non persistant for the species
    experiment={'registers':registers, 'canonicalize':canonicalize, 'tool':tool, 'nthreads':nthreads, 'debug':debug, 'baseset':set(), 'safety':safety}

    # create a SpeciesSpecifics object that will tell us where the input files can be found and keep track of where the output files should be written
    speciesinfo = SpeciesSpecifics(tag=tag, genomedir=genomedir, sketchdir=sketchdir, kstart=kstart, tool=tool, flist_loc=flist_loc)
    #inputdir = speciesinfo.inputdir
    fastas=[]
    if flist_loc:
        with open(flist_loc) as file:
            fastas = [line.strip() for line in file]
    # right now we expect that if we are provided with a genome directory AND a file list the file list will only contain the basenames
    elif os.path.exists(speciesinfo.inputdir):
        # print("I think the input directory exists")
        fastas = speciesinfo.retrieve_fasta_files(full=True)
    else:
        ValueError("You must provide either an existing directory of fastas or a file listing the paths of the desired fastas. The directory you provided was {speciesinfo.inputdir}.")
        #fastas = [f for f in allfastas if f in fastas]
    # If a fasta file list is provided subset the fastas from the species directory to only use the intersection
    # if flist_loc:
    #     with open(flist_loc) as file:
    #         fsublist = [line.strip() for line in file]
    #     fastas = [f for f in fastas if f in fsublist]
    fastas.sort()
    if nchildren:
        dtree = DeltaTree(fasta_files=fastas,speciesinfo=speciesinfo, nchildren=nchildren, experiment=experiment)
    else:
        dtree = DeltaSpider(fasta_files=fastas, speciesinfo=speciesinfo, experiment=experiment)
    # Save the cardinality keys as well as the fasta to hex dictionary lookup for the next run of the species
    
    # cardpath=os.path.join(speciesinfo.sketchdir, f'{tag}_{tool}_cardinalities.pickle')
    speciesinfo.save_cardkey(tool=tool)
    speciesinfo.save_references()
    #print(dtree)#.print_tree()
    return dtree