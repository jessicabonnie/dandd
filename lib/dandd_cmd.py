#!/usr/bin/env python
import argparse
import huffman_dandd
from huffman_dandd import write_listdict_to_csv
import sys
import pickle
import os
# from tabulate import tabulate
import csv
from typing import List, Dict, Set, Tuple, NamedTuple
# from allpairs import summ_to_phylip, rename_seqids_in_tree
# from allpairs import run_fneighbor
import json

# subcommand help: https://stackoverflow.com/questions/362426/implementing-a-command-action-parameter-style-command-line-interfaces

def insert_pre_ext(filename, string):
    toks = filename.split('.')
    return '.'.join(toks[:-1] + [string] + [toks[-1]])

def tree_command(args):
    if not (args.genomedir or args.flist_loc):
        print("ERROR: You must provide either a datadirectory or a fasta file list!")
        sys.exit(1)
    if not args.sketchdir:
        # args.sketchdir=os.path.join(args.outdir,args.tag,"sketchdb")
        args.sketchdir=os.path.join(args.outdir,"sketchdb")
        os.makedirs(args.sketchdir, exist_ok=True)
    tool="dashing"
    if args.exact:
        tool='kmc'
        args.registers=20
    if args.ksweep:
        args.ksweep=(int(args.mink), int(args.maxk))
        
    os.makedirs(args.outdir, exist_ok=True)
    # os.chdir(args.sketchdir)
    dtree = huffman_dandd.create_delta_tree(tag=args.tag, genomedir=args.genomedir, sketchdir=args.sketchdir, kstart=args.kstart, nchildren=args.nchildren, registers=args.registers, flist_loc=args.flist_loc, canonicalize=args.canonicalize, tool=tool, debug=args.debug, nthreads=args.nthreads, safety=args.safety, fast=args.fast, verbose=args.verbose, ksweep=args.ksweep, lowmem=args.lowmem)

    fileprefix = dtree.make_prefix(outdir=args.outdir, tag=args.tag, label=args.label)
    dtree.save(fileprefix=fileprefix, fast=args.fast)
   

def progressive_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
    if not args.tag:
        args.tag=dtree.speciesinfo.tag
    args.outfile = dtree.make_prefix(tag=args.tag, label=f"progu{args.norderings}", outdir=args.outdir)
        
    dtree.experiment["debug"] = args.debug
    dtree.experiment["safety"] = args.safety
    dtree.experiment["fast"] = args.fast
    dtree.experiment["verbose"] = args.verbose
    dtree.experiment["ksweep"] = None
    dtree.experiment["lowmem"] = args.lowmem
    dtree.experiment["baseset"] = set()
    if args.ksweep:
        dtree.experiment["ksweep"]=(int(args.mink), int(args.maxk))
    dtree.speciesinfo.update(tool=dtree.experiment["tool"])
    results, summary=dtree.progressive_wrapper(flist_loc=args.flist_loc, count=args.norderings, ordering_file=args.ordering_file, step=args.step)
    
    write_listdict_to_csv(outfile=args.outfile+'.csv', listdict=results)
    write_listdict_to_csv(outfile=args.outfile+'summary.csv', listdict=summary)
    # if len(str.split(args.outfile)<2):
    #     args.outfile=os.path.join(os.path.curdir,args.outfile)
    dtree.save(fileprefix=args.outfile)
    # dtree.save(outdir=os.path.split(args.outfile)[0], tag=dtree.speciesinfo.tag, label=f"progu{args.norderings}")

def abba_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
    # results=dtree.progressive_wrapper(flist_loc=args.flist_loc, count=args.norderings, ordering_file=args.ordering_file, step=args.step)
    # if args.outfile:
    #     results.to_csv(args.outfile)
    # else:
    #     print(tabulate(results, headers=list(results.columns)))

def info_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
    if not args.tag:
        args.tag=dtree.speciesinfo.tag
    if args.ksweep:
        dtree.experiment["ksweep"]=(int(args.mink), int(args.maxk))
    args.outfile = dtree.make_prefix(tag=args.tag, label=f"info", outdir=args.outdir)
    summary=dtree.summarize_tree(mink=args.mink, maxk=args.maxk)
    write_listdict_to_csv(outfile=args.outfile+'_ksweep.csv', listdict=summary)

def kij_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
    dtree.speciesinfo.update(tool=dtree.experiment["tool"])
    if not args.tag:
        args.tag=dtree.speciesinfo.tag
    args.outfile = dtree.make_prefix(tag=args.tag, label='kij', outdir=args.outdir)
    fastas=[]
    if args.flist_loc:
        with open(args.flist_loc) as file:
            fastas = [line.strip() for line in file]

    dtree.experiment["ksweep"] = (args.mink, args.maxk)
    dtree.ksweep(mink=int(args.mink),maxk=int(args.maxk))
    
    kij_results, j_results = dtree.pairwise_spiders(sublist=fastas, mink=args.mink, maxk=args.maxk, jaccard=args.jaccard)
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

    if args.jaccard:
        write_listdict_to_csv(outfile=args.outfile,listdict=j_results,suffix=".jaccard.csv")


def parse_arguments():
    parser = argparse.ArgumentParser(prog="DandD", description='program to explore delta values for a set of fasta files')
    
    # Arguments for top-level
    commands = []
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument("--verbose", "-v", action="store_true", default=False, help="Print some trees and report steps of actions.")
    parser.add_argument( "--debug", action="store_true", default=False, dest="debug", help="Share command calls to 3rd party programs.")
    parser.add_argument( "--lowmem", action="store_true", default=False, dest="lowmem", help="Delete all multi-fasta sketches, while keeping cardinality stored in dictionary for later runs. Do not recommend using with --safe")
    parser.add_argument( "--safe", action="store_true", default=False, dest="safety",   help="Double check all sketch hashes to make sure they match the sums of the fasta hashes.")
    parser.add_argument( "--fast", action="store_true", default=False, dest="fast",   help="Don't save so much stuff for second usage.")

    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands',help='additional help')
    subparsers.required = True

    # Make parser for "dand_cmd.py tree ..."
    tree_parser = subparsers.add_parser("tree", help="Calculate deltas for input fastas and full union. Create DandD tree object for further downstream analysis.")
    commands.append('tree')
    # tree_parser.add_argument( "--debug", action="store_true", default=False, dest="debug")

    tree_parser.add_argument("-s", "--tag", dest="tag", help="tagname used to label outputfiles; if datadir contains subdirectory by the same name fastas will be sourced from there",  metavar="PREFIX TAG", type=str, required=False, default='dandd')
    # choices=['ecoli', 'salmonella', 'human', 'HVSVC2','HVSVC2_snv', 'HVSVC2_snv_indel','HVSVC2_snv_sv', 'bds']
    tree_parser.add_argument("-x", "--exact", dest="exact", help="instead of estimating, count kmers using kmc3", default=False, action="store_true", required=False)

    tree_parser.add_argument("-d", "--datadir", dest="genomedir", default=None, help="data directory containing the fasta files -- all will be included if --fastas is not used", type=str, metavar="FASTADIR")

    tree_parser.add_argument("-o", "--out", dest="outdir", default=os.getcwd(), help="top level output directory that will contain the output files after running", type=str, metavar="OUTPUT DIR")

    tree_parser.add_argument("-c", "--sketchdir", dest="sketchdir", default=None, help="sketch directory to use for experiment. Default to sketchdb inside the top level output directory", type=str, metavar="SKETCHDIR")

    tree_parser.add_argument("-k", "--kstart", dest="kstart", default=12, help="kmer length at which to start the search for delta (different species have different optimal k values)", type=int, metavar="KSTART")
    
    # tree_parser.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout. NOT IMPLEMENTED.")

    tree_parser.add_argument("-f", "--fastas", dest="flist_loc", metavar="FILEPATH", type=str, default=None, help="filepath to a subset of fasta files to use in the species directory -- no title, one per line")

    tree_parser.add_argument("-l", "--label", dest="label", metavar="SUFFIX TAG", default="", help="NOT IMPLEMENTED. label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list).", required=False)

    tree_parser.add_argument("-n", "--nchildren", dest="nchildren", metavar="INTEGER", type=int, default=None, help="number of children for each node in the delta tree -- default is to create a tree of only 2 levels with all individual sketches as the children of the root node.")

    tree_parser.add_argument("-r", "--registers", dest="registers", metavar="INTEGER", default=20, help="number of registers to use during sketching")

    tree_parser.add_argument("-e", "--nthreads", dest="nthreads", metavar="INTEGER", type=int, default=10, help="number of threads to use in calls to KMC ONLY. Dashing is not currently set to use threads.")


    tree_parser.add_argument("-C", "--no-canon", action="store_false", default=True,  dest="canonicalize", help="instruct dashing to use non-canonicalized kmers")

    tree_parser.add_argument( "--ksweep", dest="ksweep", default=None,
    action="store_true", help="indicate whether a k sweep should be performed for the combinations. Without --mink and --maxk, will default to mink=2, maxk=32")
    tree_parser.add_argument("--mink", dest="mink", metavar="MINIMUM-K", required=False, default=2, type=int, help="Minimum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    tree_parser.add_argument("--maxk", dest="maxk", metavar="MAXIMUM-K", required=False, default=32, type=int, help="Maximum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    tree_parser.set_defaults(func=tree_command)


    # Make parser for "dand_cmd.py progressive ..."
    progressive_parser = subparsers.add_parser("progressive", help="Measure Delta as each individual fasta is added to the set. If a specific ordering is not provided, a set of random orderings can be generated. NOTE: Options used during creation of delta tree will be used (e.g. exact/estimate, genome directory, species tag name.)")
    commands.append('progressive')

    progressive_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="DELTA TREE", required=True, help="filepath to a pickle produced by the tree command")
    progressive_parser.add_argument("-s", "--tag", dest="tag", help="tagname used to label outputfiles, default to original tag used to create input tree",  metavar="species/experiment-tag-string", type=str, required=False)

    progressive_parser.add_argument("-r", "--orderings", dest="ordering_file", metavar="ORDERING PICKLE", type=str, default=None, help="filepath to a pickle of orderings if different from default named using tag")

    progressive_parser.add_argument("-f", "--fastas", dest="flist_loc", default=None, type=str, metavar="FILEPATH", help="filepath to a subset of fasta files from the original tree which should be analyzed using progressive union. When count is not provided, the ordering in the file will be used for a single progression. The ordering will not be added to the ordering pickle.")

    progressive_parser.add_argument("-n", "--norderings", dest="norderings", default=0, type=int, help="number of random orderings to explore. If not provided, the orderings stored in the ordering pickle will be used. If that file does not exist / is not provided, program will terminate.", metavar="NUM")

    progressive_parser.add_argument("-o", "--outdir", dest="outdir", default=os.getcwd(), type=str, help="directory to write the output tables and tree.", metavar="OUTPUT DIR")

    progressive_parser.add_argument("-l", "--label", dest="label", metavar="SUFFIX TAG", default="", help="NOT IMPLEMENTED label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list).", required=False)

    progressive_parser.add_argument( "--ksweep", dest="ksweep", default=None,
    action="store_true", help="indicate whether a k sweep should be performed for the combinations. Without --mink and --maxk, will default to mink=2, maxk=32")

    progressive_parser.add_argument("--mink", dest="mink", metavar="MINIMUM-K", required=False, default=2, type=int,help="Minimum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    progressive_parser.add_argument("--maxk", dest="maxk", metavar="MAXIMUM-K", required=False, default=32, type=int, help="Maximum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    # progressive_parser.add_argument("-o", "--out", dest="outdir", default=os.getcwd(), help="top level output directory that will contain the output files after running", type=str, metavar="OUTDIRPATH")

    progressive_parser.add_argument("--step", dest="step", default=1, type=int, help="Number of sketches to include in each progression. Mostly used for a single ordered progression.", metavar="INTEGER")

    # progressive_parser.add_argument("-e", "--exhaustive", action="store_true", dest="exhaustive", default=False, help="instead of a randome ordering, all possible permutations of the indicated fastas will be run. Warning, this should be restricted to a small subset unless space and time considerations are managed.")

    progressive_parser.set_defaults(func=progressive_command)

    # Make parser for "dand_cmd.py info ..."
    info_parser = subparsers.add_parser("info")
    commands.append('info')
    
    info_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="DELTA TREE", required=True, help="filepath to a pickle produced by the tree command. Tree nodes will be updated to hold additional sketches as needed to perform info commands selected.")
    info_parser.add_argument("-s", "--tag", dest="tag", help="tagname used to label outputfiles, default to original tag used to create input tree",  metavar="PREFIX TAG", type=str, required=False)
    
    info_parser.add_argument("--mink", dest="mink", metavar="MINIMUM-K", required=False, default=2, type=int, help="Minimum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    info_parser.add_argument("--maxk", dest="maxk", metavar="MAXIMUM-K", required=False, default=32, type=int, help="Maximum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    # info_parser.add_argument("-o", "--outfile", dest="outfile", default=None, type=str, help="path to write the output table. If path not provided, default will use the tag for the provided tree.")
    info_parser.add_argument("-o", "--outdir", dest="outdir", default=os.getcwd(), type=str, help="directory to write the output tables.", metavar="OUTPUT DIR")

    info_parser.add_argument("-l", "--label", dest="label", default="", help="NOT IMPLEMENTED Label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list).", required=False, metavar="SUFFIX TAG")

    info_parser.add_argument( "--ksweep", dest="ksweep", default=None,
    action="store_true", help="indicate whether a k sweep should be performed for the combinations. Without --mink and --maxk, will default to mink=2, maxk=32")


    info_parser.set_defaults(func=info_command)
    
    # Make parser for "dand_cmd.py abba ..."
    abba_parser = subparsers.add_parser("abba", help='NOT IMPLEMENTED. "A before B, B before A" runs runs all permutations of orderings of subsets where both fasta A and fasta B are present -- with analysis comparing those where (1) fasta A preceeds B in the ordering; (2) fasta B preceeds A in the ordering. NOTE: A and B should both be present in the provided tree.')
    commands.append('abba')

    abba_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="DELTA TREE", required=True, help="filepath to a pickle produced by the tree command")

    abba_parser.add_argument("-A", "--fastaA", dest="fastaA", metavar="FASTA_A_PATH", type=str, default=None, help="filepath of fasta A")

    abba_parser.add_argument("-B", "--fastaB", dest="fastaB", metavar="FASTA_B_PATH", type=str, default=None, help="filepath of fasta B")

    abba_parser.set_defaults(func=abba_command)

   # Make parser for "dand_cmd.py kij ..."
    kij_parser = subparsers.add_parser("kij", help="K Independent Jaccard. If a subset of fastas is not provided, matrix will include all inputs used to generate the delta tree using the `tree` command. NOTE: Options used during creation of delta tree will be used (e.g. exact/estimate, genome directory, species tag name.)")
    commands.append('kij')

    kij_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="DELTA TREE", required=True, help="filepath to a pickle produced by the tree command")
    kij_parser.add_argument("-s", "--tag", dest="tag", help="tagname used to label outputfiles, default to original tag used to create input tree",  metavar="PREFIX TAG", type=str, required=False)

    kij_parser.add_argument("-f", "--fastas", dest="flist_loc", default=None, type=str, metavar="FILEPATH", help="filepath to a subset of fasta files from the original tree which should be analyzed.")

    kij_parser.add_argument("-o", "--outdir", dest="outdir", default=os.getcwd(), type=str, help="directory to write the output tables.", metavar="OUTPUT DIR")

    kij_parser.add_argument("-l", "--label", dest="label", metavar="SUFFIX TAG", default="", help="NOT IMPLEMENTED Label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list).", required=False)

    kij_parser.add_argument("--afproject", dest="afproject", default=False, action="store_true", help="Indicate whether intermediate pickle of tuples should be created for use with helpers/afproject.py.")

    # kij_parser.add_argument("-o", "--out", dest="outdir", default=os.getcwd(), help="top level output directory that will contain the output files after running", type=str, metavar="OUTDIRPATH")

    kij_parser.add_argument("--jaccard", dest="jaccard", default=False, action="store_true", help="Indicate whether to include output for standard jaccard difference for indicated ks")

    # kij_parser.add_argument("--mash", dest="mash", default=False, action="store_true", help="Indicate whether to include output for mash difference for indicated ks")

    kij_parser.add_argument("--mink", dest="mink", metavar="MINIMUM-K", required=False, type=int, default=2, help="Minimum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    kij_parser.add_argument("--maxk", dest="maxk", metavar="MAXIMUM-K", required=False, type=int, default=32, help="Maximum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    kij_parser.set_defaults(func=kij_command)

    return parser, commands



# def main():
#     parser, commands = parse_arguments()
#     args = parser.parse_args(sys.argv[1:])
#     if len(sys.argv) < 2:
#         print('Must specify a command: ' + str(commands), file=sys.stderr)
#         return 1
#     args.func(args)
#     return 0



# if __name__ == '__main__':
#     sys.exit(main())