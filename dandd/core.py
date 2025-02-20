from __future__ import annotations
import os
import sys
import pickle
from dandd.utils import write_listdict_to_csv
# from .sketch_classes import SketchObj, DashSketchObj, KMCSketchObj, SketchFilePath
from .sketch_base import SketchObj
from .sketch_dashing import DashSketchObj
from .sketch_kmc import KMCSketchObj
from .sketch_filepath import SketchFilePath
from .species_specifics import SpeciesSpecifics
from .delta_tree import DeltaTree, DeltaSpider


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

class DandD:
    """Main class implementing DandD functionality"""
    
    def __init__(self, debug=False, fast=False, safe=False, verbose=False):
        self.debug = debug
        self.fast = fast
        self.safe = safe
        self.verbose = verbose

    def run_tree(self, args):
        """Implement tree command functionality"""
        if self.verbose:
            print("Running tree command...")
        
        if not (args.genomedir or args.flist_loc):
            print("ERROR: You must provide either a datadirectory or a fasta file list!")
            sys.exit(1)
        if not args.sketchdir:
            args.sketchdir = os.path.join(args.outdir, "sketchdb")
            os.makedirs(args.sketchdir, exist_ok=True)
        
        tool = "dashing"
        if args.exact:
            tool = 'kmc'
            args.registers = 20
        if args.ksweep:
            args.ksweep = (int(args.mink), int(args.maxk))
        
        os.makedirs(args.outdir, exist_ok=True)
        
        # Print debug info if verbose
        if self.verbose:
            print(f"Using tool: {tool}")
            print(f"Output directory: {args.outdir}")
            print(f"Sketch directory: {args.sketchdir}")
            if args.genomedir:
                print(f"Genome directory: {args.genomedir}")
            if args.flist_loc:
                print(f"Fasta list: {args.flist_loc}")
        
        # Create the delta tree with explicit parameters
        dtree = create_delta_tree(
            tag=args.tag,
            genomedir=args.genomedir,
            sketchdir=args.sketchdir,
            kstart=args.kstart,
            nchildren=args.nchildren,
            registers=int(args.registers),  # Ensure registers is an int
            flist_loc=args.flist_loc,
            canonicalize=args.canonicalize,
            tool=tool,
            debug=self.debug,
            nthreads=int(args.nthreads),
            safety=self.safe,
            fast=self.fast,
            verbose=self.verbose,
            ksweep=args.ksweep,
            lowmem=args.lowmem
        )

        if self.verbose:
            print("Delta tree created, saving results...")

        fileprefix = dtree.make_prefix(outdir=args.outdir, tag=args.tag, label=args.label)
        dtree.save(fileprefix=fileprefix, fast=self.fast)
        
        if self.verbose:
            print(f"Results saved to {fileprefix}")

    def run_progressive(self, args):
        """Implement progressive command functionality"""
        if self.verbose:
            print("Running progressive command...")
        
        dtree = pickle.load(open(args.delta_tree, "rb"))
        if not args.tag:
            args.tag = dtree.speciesinfo.tag
        args.outfile = dtree.make_prefix(tag=args.tag, label=f"progu{args.norderings}", outdir=args.outdir)
        
        dtree.experiment["debug"] = self.debug
        dtree.experiment["safety"] = self.safe
        dtree.experiment["fast"] = self.fast
        dtree.experiment["verbose"] = self.verbose
        dtree.experiment["lowmem"] = args.lowmem
        dtree.experiment["baseset"] = set()
        dtree.experiment["ksweep"] = None
        
        if args.ksweep:
            dtree.experiment["ksweep"] = (int(args.mink), int(args.maxk))
        dtree.speciesinfo.update(tool=dtree.experiment["tool"])
        
        results, summary = dtree.progressive_wrapper(
            flist_loc=args.flist_loc,
            count=args.norderings,
            ordering_file=args.ordering_file,
            step=args.step
        )
        
        write_listdict_to_csv(outfile=args.outfile + '.csv', listdict=results)
        write_listdict_to_csv(outfile=args.outfile + 'summary.csv', listdict=summary)
        dtree.save(fileprefix=args.outfile)

    def run_kij(self, args):
        """Implement k-independent Jaccard functionality"""
        if self.verbose:
            print("Running KIJ command...")

        dtree = pickle.load(open(args.delta_tree, "rb"))
        dtree.speciesinfo.update(tool=dtree.experiment["tool"])
        if not args.tag:
            args.tag=dtree.speciesinfo.tag
        args.outfile = dtree.make_prefix(tag=args.tag, label=args.label, 
        outdir=args.outdir)
        fastas=[]
        if args.flist_loc:
            with open(args.flist_loc) as file:
                fastas = [line.strip() for line in file]

        if args.ksweep:
            dtree.experiment["ksweep"]=(int(args.mink), int(args.maxk))
        dtree.ksweep(mink=int(args.mink),maxk=int(args.maxk))
        
        kij_results, j_results = dtree.pairwise_spiders(sublist=fastas, mink=args.
        mink, maxk=args.maxk, jaccard=args.jaccard)
        write_listdict_to_csv(outfile=args.outfile+".kij.csv", listdict=kij_results)
        if args.jaccard:
            write_listdict_to_csv(outfile=args.outfile+".j.csv", listdict=j_results)
        #print(j_and_kij_summ)
        dtree.speciesinfo.save_cardkey(dtree.experiment["tool"])
        dtree.speciesinfo.save_references(fast=False)
        if args.afproject:
            j_and_kij_summ = dtree.prepare_AFproject(kij_results, j_results)
            with open(args.outfile+"_AFtuples.pickle","wb") as f:
                pickle.dump(obj=j_and_kij_summ, file=f)

    # @staticmethod
    # # def create_delta_tree(tag: str, genomedir: str, sketchdir: str, kstart: int, nchildren=None, registers=0, flist_loc=None, canonicalize=True, tool='dashing', debug=False, nthreads=0, safety=False, fast=False, verbose=False, ksweep=None, lowmem=False):
    #     '''Given a species tag and a starting k value retrieve a list of fasta files to create a tree with the single fasta sketches populating the leaf nodes and the higher level nodes populated by unions
    #     tag = species tag
    #     genomedir = parent directory of species subdirectory
    #     sketchdir = parent directory where output sketches should be created
    #     kstart = starting k to use while searching for delta
    #     nchildren = number of children that nodes should have (until they can't)
    #     registers = number of registers to use when sketching
    #     flist_loc = file containing list of subset of fasta files to use from species directory (IN FUTURE maybe list of fastas with loc?)
    #     canonicalize = T/F indicating whether kmers should be canonicalized
    #     tool = string indicating which tool to use for kmer cardinality
    #     choices=["dashing","kmc"] '''
    #     # create an experiment dictionary for values that are needed at multiple levels that are non persistant for the species
    #     experiment={'registers':registers, 'canonicalize':canonicalize, 'tool':tool, 'nthreads':int(nthreads), 'debug':debug, 'baseset':set(), 'safety':safety, 'fast':fast, 'verbose':verbose, 'ksweep':ksweep, 'lowmem': lowmem}

    #     # create a SpeciesSpecifics object that will tell us where the input files can be found and keep track of where the output files should be written
    #     speciesinfo = SpeciesSpecifics(tag=tag, genomedir=genomedir, sketchdir=sketchdir, kstart=kstart, tool=tool, flist_loc=flist_loc)
        
    #     #inputdir = speciesinfo.inputdir
    #     fastas=[]
    #     if flist_loc:
    #         with open(flist_loc) as file:
    #             fastas = [line.strip() for line in file]
    #     # right now we expect that if we are provided with a genome directory AND a file list the file list will only contain the basenames
    #     elif os.path.exists(speciesinfo.inputdir):
    #         # print("I think the input directory exists")
    #         fastas = speciesinfo.retrieve_fasta_files(full=True)
    #     else:
    #         ValueError("You must provide either an existing directory of fastas or a file listing the paths of the desired fastas. The directory you provided was {speciesinfo.inputdir}.")
    #     fastas.sort()
    #     if nchildren:
    #         dtree = DeltaTree(fasta_files=fastas,speciesinfo=speciesinfo, nchildren=nchildren, experiment=experiment)
    #     else:
    #         dtree = DeltaSpider(fasta_files=fastas, speciesinfo=speciesinfo, experiment=experiment)

    #     # Save the cardinality keys as well as the fasta to hex dictionary lookup for the next run of the species
    #     speciesinfo.save_cardkey(tool=tool,fast=fast)
    #     speciesinfo.save_references(fast=fast)
    #     return dtree