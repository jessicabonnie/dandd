import sys
import os
import pickle
import csv
from random import sample, shuffle
from species_specifics import SpeciesSpecifics
from sketch_classes import *
from typing import List, Dict, Set, Tuple
from math import factorial
from itertools import permutations


def write_listdict_to_csv(outfile: str, listdict:List[Dict], suffix:str="", last_col: str = None):
    '''
    Take a list of dictionaries and write them to a csv files with the keys as headers. If a particular column needs to be last, provide the name.
    '''
    writer = open(outfile+suffix, "w") if outfile is not None and outfile != '-' else sys.stdout
    fieldnames=set()
    for x in listdict:
        fieldnames.update(x.keys())
    fieldnames=list(fieldnames)
    # The fastas field should be at the end since it sometimes has commas
    if "fastas" in fieldnames:
        last_col="fastas"
    if "files" in fieldnames:
        last_col="files"
    if last_col:
        i=fieldnames.index(last_col)
        fieldnames=fieldnames[:i]+fieldnames[i+1:]+ [fieldnames[i]]
    dict_writer = csv.DictWriter(writer, fieldnames=fieldnames)
    dict_writer.writeheader()
    dict_writer.writerows(listdict)
    writer.close()


def permute(length, norder, preexist=set(), exhaust=False, verbose=False)-> Set[Tuple[int]]:
    '''Create a set of ordering tuples. One approach will be used if the maximum number of possible permutations is low or desired. Another will be used if not.'''
    fact=factorial(length)
    newset=preexist.copy()
    norder = min(norder, fact)
    if verbose:
        print(f"{norder} permutations will be produced.")
    # if the number of permutations is low enough, generate all of them
    if fact < 5041 or norder==fact or exhaust==True:
        if verbose:
            print(f"Factorial of {norder} less than 5041 (7! + 1). Randomized list of all permutations will be produced and then subset.")
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
        #     self.ksketches = self.ksketches + [None]* (maxk-len(self.ksketches))
            self.ksketches.extend([None]* (kval-len(self.ksketches)+2))
        # make sure that all necessary ingredient sketches are available for the current k in question
        # NOTE make sure the double check which of these is necessary
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
        # seek it after a small jump also
        #self.find_delta_helper(kval = self.bestk + 2, direction=1)
        #self.find_delta_helper(kval = self.bestk + 4, direction=1)
        #TODO maybe this one should reference the species kstart
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
            # self.ksketches = self.ksketches + [None]* (maxk-len(self.ksketches))
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
            #TODO: Fiddling was happening here
            for i in range(len(self.children)):
                sketchlist = sketchlist + self.children[i].ksweep_update_node(mink=mink, maxk=maxk)
                presketches= presketches + [self.children[i].ksketches[0].sfp.full]
        # store a "sketchobj" at k=0 that holds the parallel command 
        if self.experiment["tool"] == "dashing":
            self.ksketches[0] = DashSketchObj(kval = 0, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)
        elif self.experiment["tool"] == "kmc":
            # TODO: create tmp directory for this obj to use here?
            self.ksketches[0] = KMCSketchObj(kval = 0, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)
            tmpdir=tempfile.mkdtemp()
            for k in empty_ks:
                os.mkdir(os.path.join(tmpdir,"k"+str(k)))
        # make a list of ks that don't have cardinalities/etc. stored in the tree and then figure out what the paths to those sketches would be
        for k in static_empty_ks:
            sketch_loc=self.ksketches[0].sfp.full.replace("{}",str(k))
            update_sketch = False
            # first check if those sketches exist
            # NOTE THIS IS A PROBLEM FOR LOW MEM
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
                # self.ksketches[k].card = float(self.speciesinfo.cardkey[sketch_loc])
                # self.ksketches[k].delta_pos = float(self.ksketches[k].card)/k
                # self.ksketches[k].card = float(self.speciesinfo.cardkey[sketch_loc])
                # self.ksketches[k].delta_pos =  float(self.speciesinfo.cardkey[sketch_loc])/k
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
        
        #self.speciesinfo.save_references()
        # for k in static_empty_ks:
        #     self.update_node(kval=int(k))
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
                # TODO can be done in parallel
                for i in range(len(self.children)):
                    self.children[i].update_node(kval)
                    presketches= presketches + [self.children[i].ksketches[kval].sketch]
            # create a sketch object dependent on whether the tool is dashing or kmc
            if self.experiment["tool"] == "dashing":
                self.ksketches[kval] = DashSketchObj(kval = kval, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)
            elif self.experiment["tool"] == "kmc":
                # TODO: create tmp directory for this obj to use here?
                self.ksketches[kval] = KMCSketchObj(kval = kval, sfp = sfp, speciesinfo=self.speciesinfo, experiment=self.experiment, presketches=presketches)
        return

    def node_ksweep(self, mink, maxk):
        '''Sketch all of the ks for the node (and its decendent nodes)between mink and maxk (even when they weren't needed to calculate delta'''
        sketchlist=self.ksweep_update_node(mink=mink, maxk=maxk)
        multicard=self.ksketches[0].card_command(sketchlist)
        self.ksketches[0].individual_card(cmd=multicard)

        # This loop is necessary to populate properties of the sketch objects
        for kval in range(mink, maxk+1):
            if self.ksketches[kval] is None:
                self.update_node(kval)
            # if self.ksketches[kval].card == 0:
            #     print("SKETCH CARD NOT THERE?")
                # if not self.ksketches[kval]:
                #     exit(f"Something is wrong with this node's k={kval} sketch. Recommend deleting the input sketches: ", self.node_title)
        self.mink=mink
        self.maxk=maxk
        return 

    def summarize(self, mink:int=0, maxk:int=0, ordering_number=0):
        '''create a dataframe of all "possible" delta values that were examined during creation of the node for use in plotting -- this isn't actually summarizing, so the function should be renamed'''
        nodevals=[]
        
        # self.node_ksweep(mink=mink,maxk=maxk)
        for kval in range(mink, maxk+1):
        #     if self.ksketches[kval] is None or self.ksketches[kval].delta_pos == 0:
        #         self.update_node(kval)
            linedict = {"ngen": self.ngen, "kval": kval, "card": self.ksketches[kval].card, "delta_pos": self.ksketches[kval].delta_pos, "title": self.node_title, "command": self.ksketches[kval].cmd, "ordering":ordering_number}
            # for index, value in enumerate(self.fastas):
            #     linedict[f"step{index}"] = value
            nodevals.append(linedict)
        return nodevals


class DeltaTree:
    ''' Delta tree data structure. '''
    def __init__(self, fasta_files, speciesinfo, nchildren=2, leafnodes=[], experiment={'tool':'dashing', 'registers':20, 'canonicalize':True, 'debug':False, 'nthreads':0, 'baseset': set(), 'safety': False, 'fast': False, 'verbose': False, 'ksweep': None, 'lowmem': False}, padding=True):
        self.experiment=experiment
        # self._symbols = []
        self.mink=0
        self.maxk=0
        if self.experiment["ksweep"] is not None:
            (self.mink, self.maxk) = self.experiment["ksweep"]
        
        
        self.kstart=speciesinfo.kstart
        self.speciesinfo=speciesinfo
        if self.experiment["verbose"]:
            print("Now making tree for fastas: " + ", ".join(fasta_files))
        self._build_tree(fasta_files, nchildren)
        self.fill_tree(padding=padding)
        self.ngen = len(fasta_files)
        self.root=self._dt[-1]
        self.delta = self.root_delta()
        self.fastas = fasta_files
        if self.experiment["ksweep"] is None:
            speciesinfo.kstart = self.root_k()
        self.speciesinfo.save_references(fast=experiment['fast'])
        self.speciesinfo.save_cardkey(tool=self.experiment["tool"])
        # self.delete_sketches()
    def __sub__(self, other):
        # sub = subtraction
        print("Larger Tree Delta: ", self.delta)
        print("Subtree Delta: ", other.delta)
        print("Subtraction Result: ", self.delta - other.delta)
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

    
    def root_delta(self):
        '''retrieve the root node of the tree'''
        root=self._dt[-1]
        return root.delta
    def root_k(self):
        '''retrieve the argmax k for the root node'''
        root=self._dt[-1]
        return root.bestk

    def delete_sketches(self):
        for node in self._dt[:-1]:
            sketches=[i for i in node.ksketches if i is not None]
            if node.ngen > 1:
                for sketch in sketches:
                    sketch.remove_sketch()


    def _build_tree(self, symbol: list, nchildren: int, leafnodes: List[DeltaTreeNode] = []) -> None:
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
        if len(leafnodes) == 0:
            inputs = [
                DeltaTreeNode(
                    node_title=s, children=[], speciesinfo=self.speciesinfo, experiment=self.experiment, progeny=[]
                ) for s in symbol]
        else:
            inputs = leafnodes
        inputs.sort()
        for n in inputs:
            if self.experiment["ksweep"] is None:
                n.find_delta(self.speciesinfo.kstart)
            else:
                # NOTE Maybe make sure delta_pos is assigned here
                n.node_ksweep(mink=self.mink, maxk=self.maxk)
        
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
            # NOTE check to make sure these are necessary here
            if self.experiment["ksweep"] is None:
                new_node.find_delta(self.speciesinfo.kstart)
            else:
                new_node.node_ksweep(mink=self.mink, maxk=self.maxk)

            while idx_insert < len(self._dt)-increment and self._dt[idx_insert+increment].ngen <= new_node.ngen:
                idx_insert += increment
            
            self._dt = self._dt[:idx_insert+increment] + [new_node] + self._dt[idx_insert+increment:]
                
            idx_current += nchildren
            if idx_insert + increment > len(self._dt)-1:
                nchildren = len(self._dt) - idx_current


    def print_list(self) -> None:
        nodes = []
        for i, node in enumerate(self._dt):
            nodes.append(f'\'{node.node_title}\'({node.ngen}\'({" ".join([i.node_title for i in node.progeny])})')
        print(' -> '.join(nodes))

    def fill_tree(self, padding=False) -> None:
        '''Starting at the root make sure that all nodes in the tree contain the sketches for the argmax ks for every node as well as 2 less than the minimum and 2 greater than the maximum (IF padding argument is True)'''
        root = self._dt[-1]
        if self.experiment["ksweep"] is None:
            bestks = list(set([n.bestk for n in self._dt]))
            bestks = [k for k in bestks if k!=0  ]
            bestks.sort()
        # if padding:
        #     bestks = bestks + [bestks[0]-1] + [bestks[0]-2] + [bestks[-1]+1] + [bestks[-1]+2] + [bestks[-1]+3]
            for k in bestks:
                root.update_node(k)
        else :
            self.ksweep(mink=self.experiment["ksweep"][0], maxk=self.experiment["ksweep"][1])
        return
        #self.speciesinfo.save_references()
        #self.speciesinfo.save_cardkey(tool=self.experiment["tool"])
    
    def leaf_nodes(self) -> List[DeltaTreeNode]:
        return [child for child in self._dt if child.ngen==1]

    # def to_spider(self):
    #     '''Transform tree into a spider if it isn't'''

    #     children=self.leaf_nodes()
    #     progeny=[p.progeny for p in children]
    #     #flatten the progeny list
    #     progeny=[item for sublist in progeny for item in sublist]
    #     child_titles=[c.node_title for c in children]
    #     body_node = DeltaTreeNode(
    #         node_title="_".join(child_titles), speciesinfo=self.speciesinfo,
    #         children = children,
    #         progeny=progeny,
    #         experiment=self.experiment
    #         )
    #     body_node.find_delta(kval=self.speciesinfo.kstart)
    #     self._dt = children + [body_node]
    
    def make_prefix(self, tag: str, label="", outdir:str=None):
        if not outdir:
            outdir=os.getcwd()
        if not label == "":
            label = "_"+label
        fileprefix=os.path.join(outdir,  tag + label + "_" + str(self.ngen) + "_" + self.experiment["tool"] )
        return fileprefix
    def save(self, fileprefix:str, fast=False):
        '''
        Save the delta tree for future retrieval
        '''
        
        filepath=fileprefix + '_dtree.pickle'
        if not fast:
            with open(filepath,"wb") as f:
                pickle.dump(obj=self, file=f)
            print("Tree Pickle saved to: "+filepath)

            expmaploc=fileprefix + '_sketchdb.txt'
            explist=[self.speciesinfo.sketchinfo[item] for item in list(self.experiment["baseset"])]
            write_listdict_to_csv(outfile=expmaploc, listdict=explist)
            print(f"Output Sketch/DB mapping saved to {expmaploc}.")
        deltapath=fileprefix + '_deltas.csv'
        write_listdict_to_csv(deltapath,self.report_deltas())
        print("Deltas saved to: " + deltapath)
        return filepath

    def report_deltas(self) -> List[dict]:
        ''' Traverse the DeltaTree to return a dataframe with the delta values of the nodes in the tree.'''
        root = self._dt[-1]
        def _delta_recursive(node) -> List[Dict]:
            tmplist=[{"delta": node.delta, "k": node.bestk, "title": node.node_title, "ngen": node.ngen,"sketchloc": node.ksketches[node.bestk].sketch, "card":node.ksketches[node.bestk].card ,"fastas": "|".join(node.fastas)}]
            if node.children:
                nchild=len(node.children)
                for i in range(nchild):
                    n=node.children[i]
                    tmplist.extend(_delta_recursive(n))
            return tmplist
        
        dictlist= _delta_recursive(root)
        return dictlist
  
    # def summarize_tree(self, mink=0, maxk=0) -> List[dict]:
    #     ''' Traverse the DeltaTree to return a dataframe with all possible delta values. -- this isn't actually summarizing, so the function should be renamed'''
    #     root = self._dt[-1]
    #     if mink == 0 or maxk == 0:
    #         mink, maxk = self.mink, self.maxk #experiment["ksweep"]
    #         # mink=self.mink
    #         # if self.mink > 4:
    #         #     mink=self.mink - 2
    #     # if maxk == 0:
    #     #     maxk = self.maxk
    #         # if self.maxk <= 30:
    #         #     maxk=self.maxk + 2
    #     self.ksweep(mink=mink, maxk=maxk)
    #     sum_listdict=[]
    #     sum_listdict.extend(root.summarize(mink=mink, maxk=maxk))
    #     def _delta_pos_recursive(node) -> List[dict]:
    #         tmplist=node.summarize(mink=mink, maxk=maxk)
    #         if node.children:
    #             nchild=len(node.children)
    #             if nchild == self.ngen and self.ngen > 1:
    #                 n=node.children[self.ngen-1]
    #                 tmplist.extend(_delta_pos_recursive(n))
    #         return sum_listdict
        
    #     dictlist= _delta_pos_recursive(root)
    #     return dictlist

    def nodes_from_fastas(self, fasta_list):
        '''
        Provided a list of fastas retrieve the leaf nodes formed from those fastas
        '''
        return [node for node in self.leaf_nodes() if node.fastas[0] in fasta_list]

    def find_delta_delta(self, fasta_subset: List[str]) -> float:
        '''Provided a list of fastas in a subset, find the delta-delta values between the whole spider and a spider without the provided fastas'''
        # create list of fastas that are in the original spider that are not in the subset provided --> i.e. the complement
        fastas = [f for f in self.fastas if f not in fasta_subset]
        small_spider = SubSpider(leafnodes=self.nodes_from_fastas(fastas), speciesinfo=self.speciesinfo, experiment=self.experiment)
        print("Full Tree Delta: ", self.delta)
        print("Subtree Delta: ", small_spider.delta)
        return self - small_spider
   
    def ksweep(self, mink, maxk) -> None:
        for node in self._dt:
            node.node_ksweep(mink=mink, maxk=maxk)
        # self.speciesinfo.save_references()
        #self.speciesinfo.save_cardkey(tool=self.experiment["tool"])

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

        # if there is no ordering file at the location, time to make one    
        else:
            # if count is not provided, how will we know how many to make??
            if count<1:
                raise ValueError("You must provide a value for count when there is no default ordering file")
        
        orderings = permute(length=len(fastas), norder=count, preexist=orderings, verbose=self.experiment['verbose'])
        # save the orderings for use next run of species 
        with open(ordering_file,"wb") as f:
            pickle.dump(orderings, f)
        return fastas, list(orderings)

    def progressive_wrapper(self, flist_loc=None, count=30, ordering_file=None,step=1, debug=False)-> List[dict]:
        fastas, orderings = self.orderings_list( ordering_file=ordering_file, flist_loc=flist_loc, count=count)

        return self.progressive_union(flist=fastas, orderings=orderings, step=step)

    def progressive_union(self, flist, orderings, step) -> Tuple[List[dict], List[dict]]:
        '''Create (or use if provided) a series of random orderings to use when adding the individual fasta sketches to a union. Outputs a table with the delta values and associated ks at each stage'''

        # create a sketch of the full union of the fastas
        smain = DeltaSpider(fasta_files=flist, speciesinfo=self.speciesinfo, experiment=self.experiment)
        results=[]
        summary=[]
        for i in range(0,len(orderings)):
            if self.experiment["verbose"]:
                print(f"Now sweeping for ordering {i+1}")
            oresults, osummary = [], []
            oresults, osummary = smain.sketch_ordering(orderings[i], ordering_number=i+1, step=step)
            # for o in osummary:
            #     o["ordering"] = i+1
            results.extend(oresults)
            summary.extend(osummary)
            self.speciesinfo.save_references(fast=self.experiment['fast'])
            self.speciesinfo.save_cardkey(tool=self.experiment["tool"],fast=self.experiment['fast'])
        return results, summary

    def sketch_ordering(self, ordering, ordering_number, step=1) -> Tuple[List[dict], List[dict]]:
        '''Provided an ordering for the fastas in a tree, create sketches of the subsets within that ordering and report the deltas'''
        flen=len(ordering)
        output=[]
        summary=[]
        krange = self.experiment["ksweep"]
        if krange is None:
            krange=(self.mink, self.maxk)
        
        for i in range(1,flen+1):
            if i % step == 0:
                sublist=[self.fastas[j] for j in ordering[:i]]
                ospider=SubSpider(leafnodes=self.nodes_from_fastas(sublist), speciesinfo=self.speciesinfo, experiment=self.experiment)
                # NOTE why do we need ksweep here if we just did it during fill_tree? or did we. Well, it doesn't work without it so.
                ospider.ksweep(mink=int(krange[0]), maxk=int(krange[1]))
                output.append({"ngen":i, "kval":ospider.root_k(), "delta": ospider.delta, "ordering": ordering_number, "fastas": sublist})
                newsum = ospider.root.summarize(mink=int(krange[0]), maxk=int(krange[1]), ordering_number=ordering_number)
                summary.extend(newsum)

        return output, summary


    def pairwise_spiders(self, sublist=[], mink=0, maxk=0, jaccard=True) -> Tuple[List[dict], List[dict]]:
        '''Create values for K-Independent-Jaccard. (Two Legged Spiders) '''
        # super_spider=self.to_spider()

        if len(sublist)==0:
            sublist=self.leaf_nodes()
        pairings=[[a, b] for idx, a in enumerate(sublist) for b in sublist[idx + 1:]]
        kij_results=[]
        j_results=[]
        new_experiment = self.experiment.copy()
        new_experiment.update({'fast': True , 'safe': False, 'ksweep':None})
        if jaccard and (mink == 0 or maxk == 0) :
            if self.experiment["ksweep"]:
                (mink, maxk) = self.experiment["ksweep"]
                print("WARNING: If EITHER minimum OR maximum k are not provided with --mink and --maxk flags, DandD will default to the --ksweep values embedded in the delta-tree input.")
            else:
                print("WARNING: If BOTH minimum AND maximum k are not provided either by the input delta-tree or using --mink and --maxk, the --jaccard flag will be ignored.")
                jaccard = False
        for pair in pairings:
            pspider=SubSpider(leafnodes=pair,speciesinfo=self.speciesinfo,experiment=new_experiment)
            pspider.root.find_delta(self.root_k())
            # pspider.ksweep(mink=mink, maxk=maxk)
            kij_results.append(pspider.kij_summarize())
            if jaccard:
                pspider.ksweep(mink=mink, maxk=maxk)
                j_results.extend(pspider.jaccard_summarize(mink=mink, maxk=maxk))  
              
        # self.speciesinfo.save_references(fast=self.experiment['fast'])
        # self.speciesinfo.save_cardkey(tool=new_experiment["tool"],fast=new_experiment['fast'])
        return kij_results, j_results

    def prepare_AFproject(self, kijsummary, jsummary) -> List[Tuple]:
        '''
        Transform the kij and jaccard listdicts into the format expected by scripts in the helper folder.
        '''
        all_out=set()
        # Records in j_and_kij_summ are of the form (tool, name1, name2, k, j, k1, k2, k12)
        #  - When name1 == name2, the tuple describes a single dataset
        #  - When name1 != name2, the tuple describes a pair of datasets
        
        tool = self.experiment["tool"]
        for dictitem in kijsummary:
            outtuple_list = [
            (tool, dictitem["Atitle"], dictitem["Btitle"], 0,  dictitem["KIJ"], dictitem["Ak"], dictitem["Bk"], dictitem["ABk"])
            ]
            all_out.update(outtuple_list)

        for dictitem in jsummary:
            outtuple_list = [
                (tool, dictitem["Atitle"], dictitem["Btitle"], dictitem["kval"], dictitem["jaccard"], None, None, None)
            ]
            all_out.update(outtuple_list)
        return list(all_out)
         #  - When the record describes a J, then there are k tuples for each pair/singleton
        #    + k is positive
        #    + k1, k2, and k12 are all None
        #    + j = J (or J_k)
        #

    

class SubSpider(DeltaTree):
    def __init__(self,leafnodes,speciesinfo,experiment):
        self.speciesinfo=speciesinfo
        self.fastahex = self.speciesinfo.fastahex
        self.experiment=experiment
        # self._symbols = []
        self.kstart=self.speciesinfo.kstart
        if self.experiment["ksweep"] is not None:
            self.mink, self.maxk = self.experiment["ksweep"]
        # else:
        #     self.mink, self.maxk = self.speciesinfo.kstart, self.speciesinfo.kstart
        self._build_tree(leafnodes)
        self.root=self._dt[-1]
        self.fastas=self.root.fastas
        self.ngen = len(self.fastas)
        self.delta = None
        # NOTE : this didn't used to be here, maybe it breaks or slows?
        self.fill_tree()
        if self.experiment["ksweep"] is None:
            self.delta = self.root_delta()
        else:
            (self.mink, self.maxk) = self.experiment["ksweep"]
            # self.fill_tree()

    def _build_tree(self, leafnodes):
        children=leafnodes
        progeny=[p.progeny for p in children]
        #flatten the progeny list
        progeny=[item for sublist in progeny for item in sublist]
        child_titles=[os.path.basename(c.node_title) for c in children]
        body_node = DeltaTreeNode(
            node_title="_".join(child_titles), speciesinfo=self.speciesinfo,
            children = children,
            progeny=progeny,
            experiment=self.experiment
            )
        # NOTE: Make sure this is necessary ... maybe should be done with fill_tree
        if self.experiment["ksweep"] is None:
            body_node.find_delta(kval=self.speciesinfo.kstart)
        else:
            body_node.node_ksweep(mink=self.mink, maxk=self.maxk)
        self.mink=body_node.mink
        self.maxk=body_node.maxk
        self._dt = children + [body_node]

    def kij_summarize(self) -> Dict:
        '''Calculate k independent jaccard for the subspider'''
        ##TODO: Write this to handle more than 2 children?
        if len(self.fastas) != 2:
            raise ValueError("KIJ can only be calculated on spider/trees with 2 children")
        self.root.update_node(self.root.bestk)
        childA=self._dt[0]
        childB=self._dt[1]
        childA.update_node(childA.bestk)
        childB.update_node(childB.bestk)
        # sort in lexigraphical order so duplicates
        names = [childA.node_title,childB.node_title]
        if names != sorted(names):
            childA=self._dt[1]
            childB=self._dt[0]
        outdict={"A":childA.fastas[0], "B":childB.fastas[0],
            "Adelta":childA.delta, "Bdelta":childB.delta,
            "Ak":childA.bestk, "Bk":childB.bestk,
            "ABdelta":self.root.delta, "ABk":self.root.bestk, "Atitle": childA.node_title, "Btitle": childB.node_title}
        outdict["KIJ"]=(outdict["Adelta"] + outdict["Bdelta"]-outdict["ABdelta"])/outdict["ABdelta"]
        
        return outdict


    def jaccard_summarize(self, mink=2, maxk=32) -> List[Dict]:
        '''Calculate jaccard distance for the subspider'''
        ##TODO: Write this to handle more than 2 children?
        if len(self.fastas) != 2:
            raise ValueError("KIJ can only be calculated on spider/trees with 2 or more children")
        childA=self._dt[0]
        childB=self._dt[1]
        jevals=[]
        if [childA.node_title, childB.node_title] != [childA.node_title, childB.node_title]:
            childA = self._dt[1]
            childB = self._dt[0]
        outdict={"A":childA.fastas[0], "B":childB.fastas[0],
        "Atitle": childA.node_title, "Btitle": childB.node_title}
        self.ksweep(mink=mink,maxk=maxk)
        for k in range(mink,maxk+1):
            odict = outdict.copy()
            odict.update({"kval":k , "Acard": childA.ksketches[k].card, "Bcard": childB.ksketches[k].card, "ABcard":self.root.ksketches[k].card})
            odict["jaccard"] = (odict["Acard"] + odict["Bcard"] - odict["ABcard"])/odict["ABcard"]
            jevals.append(odict)
        return jevals


class DeltaSpider(DeltaTree):
    '''Create a structure with all single sketches in terminal nodes tied to a single union node for all of them'''
    def __init__(self, fasta_files, speciesinfo, experiment, padding=False):
        nchildren=len(fasta_files)
        super().__init__(fasta_files=fasta_files, speciesinfo=speciesinfo, experiment=experiment, nchildren=nchildren, padding=padding)
    def __init2__(self, tree:DeltaTree):
        raise NotImplementedError("initialization of spider by tree not yet implemented")
    #TODO: add a function to the superclass that adds union nodes to the _dt and then create one here that adds just the spider body. That way the function can be passed the list of childnodes and then repurposed to instantiate a spider using only child nodes

    ## TODO: init function receives a delta tree and creates a spider out of it's nodes without creating new child nodes

# class ProgressiveUnion:
#     '''Unimplemented Progressive Union Class to hold functions and objects relating to that capability'''
#     def __init__(self, deltatree: DeltaTree, orderings: list):
#         self.orderings=[]
#         self.dtree=deltatree
#         self.orderings=orderings

#         raise NotImplementedError


def create_delta_tree(tag: str, genomedir: str, sketchdir: str, kstart: int, nchildren=None, registers=0, flist_loc=None, canonicalize=True, tool='dashing', debug=False, nthreads=0, safety=False, fast=False, verbose=False, ksweep=None, lowmem=False):
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
    experiment={'registers':registers, 'canonicalize':canonicalize, 'tool':tool, 'nthreads':int(nthreads), 'debug':debug, 'baseset':set(), 'safety':safety, 'fast':fast, 'verbose':verbose, 'ksweep':ksweep, 'lowmem': lowmem}

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
    fastas.sort()
    if nchildren:
        dtree = DeltaTree(fasta_files=fastas,speciesinfo=speciesinfo, nchildren=nchildren, experiment=experiment)
    else:
        dtree = DeltaSpider(fasta_files=fastas, speciesinfo=speciesinfo, experiment=experiment)

    # Save the cardinality keys as well as the fasta to hex dictionary lookup for the next run of the species
    speciesinfo.save_cardkey(tool=tool,fast=fast)
    speciesinfo.save_references(fast=fast)
    return dtree