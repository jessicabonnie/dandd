import os
import sys
from dandd.huffman_dandd import create_delta_tree
import pickle
from dandd.utils import write_listdict_to_csv

class DandD:
    """Main class implementing DandD functionality"""
    
    def __init__(self, debug=False, fast=False, safe=False, verbose=False):
        self.debug = debug
        self.fast = fast
        self.safe = safe
        self.verbose = verbose

    def run_tree(self, args):
        """Implement tree command functionality"""
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
        
        dtree = create_delta_tree(
            tag=args.tag,
            genomedir=args.genomedir,
            sketchdir=args.sketchdir,
            kstart=args.kstart,
            nchildren=args.nchildren,
            registers=args.registers,
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

        fileprefix = dtree.make_prefix(outdir=args.outdir, tag=args.tag, label=args.label)
        dtree.save(fileprefix=fileprefix, fast=self.fast)

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