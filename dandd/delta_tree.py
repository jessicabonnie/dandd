from __future__ import annotations
import os
import pickle
from dandd.delta_node import DeltaNode
from dandd.utils import permute, write_listdict_to_csv
from typing import List, Dict, Set, Tuple


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


    def _build_tree(self, symbol: list, nchildren: int, leafnodes: List[DeltaNode] = []) -> None:
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
                DeltaNode(
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
            new_node = DeltaNode(
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
    
    def leaf_nodes(self) -> List[DeltaNode]:
        return [child for child in self._dt if child.ngen==1]

    # def to_spider(self):
    #     '''Transform tree into a spider if it isn't'''

    #     children=self.leaf_nodes()
    #     progeny=[p.progeny for p in children]
    #     #flatten the progeny list
    #     progeny=[item for sublist in progeny for item in sublist]
    #     child_titles=[c.node_title for c in children]
    #     body_node = DeltaNode(
    #         node_title="_".join(child_titles), speciesinfo=self.speciesinfo,
    #         children = children,
    #         progeny=progeny,
    #         experiment=self.experiment
    #         )
    #     body_node.find_delta(kval=self.speciesinfo.kstart)
    #     self._dt = children + [body_node]
    
    def make_prefix(self, tag: str, label="", outdir:str=None):
        '''
        Make a prefix for the output file
        '''
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

    def orderings_list(self, fastas: List[str], ordering_file=None, flist_loc=None, count=0, verbose=False)-> Tuple[List[str], List[Tuple[int]]]:
        '''create or retrieve a series of random orderings of fasta sketches. return also the expected "sorted" array of the files. A subset of the fastas in the tree can be provided by name (in a file). The ordering of this file will be used when count=1 and the list is provided.'''
        
        # if count is one "sorted" ordering is returned with the reference list
        if count == 1:
            return [tuple(i for i in range(len(fastas)))]
        
        # orderings are handled as sets to prevent duplication
        orderings=set()
        # if the ordering file exists then read the orderings
        if os.path.exists(ordering_file):
            with open(ordering_file,'rb') as f:
                orderings=pickle.load(f)
            # if count was not provided then just return the orderings that are already there
            if count==0:
                return list(orderings)
            # if the count is lte to the number of orderings in the file, take the first count number of orderings
            if count <= len(orderings):
                return list(orderings)[:count]

        # if there is no ordering file at the location, time to make one    
        else:
            # if count is not provided, how will we know how many to make??
            if count<1:
                raise ValueError("You must provide a value for count when there is no default ordering file")
        
        orderings = permute(length=len(fastas), norder=count, preexist=orderings, verbose=verbose)
        # save the orderings for use next run of species 
        with open(ordering_file,"wb") as f:
            pickle.dump(orderings, f)
        return list(orderings)

    def progressive_wrapper(self, flist_loc=None, count=30, ordering_file=None,step=1, debug=False)-> List[dict]:
        '''
        Wrapper for the progressive union function.
        '''
        fastas = self.subset_fastas(flist_loc=flist_loc)
        # if no ordering file is provided the default location is used
        if not ordering_file:
            default_ordering=os.path.join(self.speciesinfo.sketchdir, self.speciesinfo.tag + "_"+ str(len(fastas))+"_orderings.pickle")
            ordering_file=default_ordering

        orderings = self.orderings_list(fastas=fastas, ordering_file=ordering_file, flist_loc=flist_loc, count=count, verbose=self.experiment['verbose'])

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


    def subset_fastas(self, flist_loc: str) -> List[str]:
        '''
        Provided a file of fasta names, return a list of fastas that are in the species directory.
        '''
        fastas=self.fastas
        fastas.sort()

        # If a fasta file list is provided, subset the fastas from the species directory to only use the intersection
        if flist_loc:
            with open(flist_loc) as file:
                fsublist = [line.strip() for line in file]
            fastas = [f for f in fastas if f in fsublist]
            #order fastas as given in file
            fastas = [f for f in fsublist if f in fastas]
        return fastas

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
        body_node = DeltaNode(
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
