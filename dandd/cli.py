#!/usr/bin/env python3
from __future__ import annotations
import argparse
import sys
from .core import DandD
import os

def add_universal_cmds(subparser:argparse.ArgumentParser):
    subparser.add_argument('--version', action='version', version='%(prog)s 1.0.0')
    subparser.add_argument("--verbose", "-v", action="store_true", default=False, help="Print some trees and report steps of actions.")
    subparser.add_argument( "--debug", action="store_true", default=False, dest="debug", help="Share command calls to 3rd party programs.")
    subparser.add_argument( "--lowmem", action="store_true", default=False, dest="lowmem", help="Delete all multi-fasta sketches, while keeping cardinality stored in dictionary for later runs. Don't require sketch/db existence if the cardinality is already stored. Do not recommend using with --safe")
    subparser.add_argument( "--safe", action="store_true", default=False, dest="safety",   help="Double check all sketch/db name hashes to make sure they match the sums of the component fasta hashes.")
    subparser.add_argument( "--fast", action="store_true", default=False, dest="fast",   help="Don't save so much stuff for second usage.")
    return subparser



# def main_parser_command(args):
#     experiment = {}
#     experiment["debug"] = args.debug
#     experiment["safety"] = args.safety
#     experiment["fast"] = args.fast
#     experiment["verbose"] = args.verbose
   


def tree_command(args):
    dandd = DandD(
        debug=args.debug,
        fast=args.fast,
        safe=args.safety,
        verbose=args.verbose
    )
    dandd.run_tree(args)

def progressive_command(args):
    dandd = DandD(
        debug=args.debug,
        fast=args.fast,
        safe=args.safety,
        verbose=args.verbose
    )
    dandd.run_progressive(args)

def kij_command(args):
    dandd = DandD(
        debug=args.debug,
        fast=args.fast,
        safe=args.safety,
        verbose=args.verbose
    )
    dandd.run_kij(args)


def parse_arguments():

     # Arguments shared across all commands
    parent_parser=argparse.ArgumentParser(add_help=False)
    parent_parser=add_universal_cmds(parent_parser)

    # ksweep argument parser
    ksweep_parser=argparse.ArgumentParser(add_help=False)
    ksweep_parser.add_argument( "--ksweep", dest="ksweep", default=None,
    action="store_true", help="indicate whether a k sweep should be performed for the combinations. Without --mink and --maxk, will default to mink=2, maxk=32")
    ksweep_parser.add_argument("--mink", dest="mink", metavar="MINIMUM-K", required=False, default=2, type=int, help="Minimum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")
    ksweep_parser.add_argument("--maxk", dest="maxk", metavar="MAXIMUM-K", required=False, default=32, type=int, help="Maximum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")


    # Top level parser
    parser = argparse.ArgumentParser(prog="DandD", 
    description='program to explore delta values for a set of fasta files',parents=[parent_parser])
    
    # Keep track of subcommands 
    commands = []
    # Create subcommand parser
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands',help='additional help',dest='command')
    subparsers.required = True

    # Make parser for "dand_cmd.py tree ..."
    tree_parser = subparsers.add_parser("tree", help="Calculate deltas for input fastas and full union. Create DandD tree object for further downstream analysis.", parents=[parent_parser, ksweep_parser])
    commands.append('tree')

    tree_parser.add_argument("-s", "--tag", dest="tag", help="tagname used to label outputfiles; if datadir contains subdirectory by the same name fastas will be sourced from there",  metavar="PREFIX TAG", type=str, required=False, default='dandd')
 
    tree_parser.add_argument("-x", "--exact", dest="exact", help="instead of estimating, count kmers using kmc3", default=False, action="store_true", required=False)

    tree_parser.add_argument("-d", "--datadir", dest="genomedir", default=None, help="data directory containing the fasta files -- all will be included if --fastas is not used", type=str, metavar="FASTADIR")

    tree_parser.add_argument("-o", "--out", dest="outdir", default=os.getcwd(), help="top level output directory that will contain the output files after running", type=str, metavar="OUTPUT DIR")

    tree_parser.add_argument("-c", "--sketchdir", dest="sketchdir", default=None, help="sketch directory to use for experiment. Default to sketchdb inside the top level output directory", type=str, metavar="SKETCHDIR")

    tree_parser.add_argument("-k", "--kstart", dest="kstart", default=12, help="kmer length at which to start the search for delta (different species have different optimal k values)", type=int, metavar="KSTART")
    

    tree_parser.add_argument("-f", "--fastas", dest="flist_loc", metavar="FILEPATH", type=str, default=None, help="filepath to a subset of fasta files to use in the species directory -- no title, one per line")

    tree_parser.add_argument("-l", "--label", dest="label", metavar="SUFFIX TAG", default="", help="NOT IMPLEMENTED. label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list).", required=False)

    tree_parser.add_argument("-n", "--nchildren", dest="nchildren", metavar="INTEGER", type=int, default=None, help="number of children for each node in the delta tree -- default is to create a tree of only 2 levels with all individual sketches as the children of the root node.")

    tree_parser.add_argument("-r", "--registers", dest="registers", metavar="INTEGER", default=20, help="number of registers to use during sketching")

    tree_parser.add_argument("-e", "--nthreads", dest="nthreads", metavar="INTEGER", type=int, default=0, help="number of threads to use in calls to KMC ONLY. Dashing is not currently set to use threads. Default uses max available cores.")

    tree_parser.add_argument("-C", "--no-canon", action="store_false", default=True,  dest="canonicalize", help="instruct dashing to use non-canonicalized kmers")

    tree_parser.set_defaults(func=tree_command)

    # Make parser for "dand_cmd.py progressive ..."
    progressive_parser = subparsers.add_parser("progressive", help="Measure Delta as each individual fasta is added to the set. If a specific ordering is not provided, a set of random orderings can be generated. NOTE: Options used during creation of delta tree will be used (e.g. exact/estimate, genome directory, species tag name.)", 
    parents=[parent_parser, ksweep_parser])
    commands.append('progressive')

    progressive_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="DELTA TREE", required=True, help="filepath to a pickle produced by the tree command")
    progressive_parser.add_argument("-s", "--tag", dest="tag", help="tagname used to label outputfiles, default to original tag used to create input tree",  metavar="species/experiment-tag-string", type=str, required=False)

    progressive_parser.add_argument("-r", "--orderings", dest="ordering_file", metavar="ORDERING PICKLE", type=str, default=None, help="filepath to a pickle of orderings if different from default named using tag")

    progressive_parser.add_argument("-f", "--fastas", dest="flist_loc", default=None, type=str, metavar="FILEPATH", help="filepath to a subset of fasta files from the original tree which should be analyzed using progressive union. When count is not provided, the ordering in the file will be used for a single progression. The ordering will not be added to the ordering pickle.")

    progressive_parser.add_argument("-n", "--norderings", dest="norderings", default=0, type=int, help="number of random orderings to explore. If not provided, the orderings stored in the ordering pickle will be used. If that file does not exist / is not provided, program will terminate.", metavar="NUM")

    progressive_parser.add_argument("-o", "--outdir", dest="outdir", default=os.getcwd(), type=str, help="directory to write the output tables and tree.", metavar="OUTPUT DIR")

    progressive_parser.add_argument("-l", "--label", dest="label", metavar="SUFFIX TAG", default="", help="NOT IMPLEMENTED label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list).", required=False)

    progressive_parser.add_argument("--step", dest="step", default=1, type=int, help="Number of sketches to include in each progression. Mostly used for a single ordered progression.", metavar="INTEGER")

    progressive_parser.set_defaults(func=progressive_command)

    # Make parser for "dand_cmd.py info ..."
    # info_parser = subparsers.add_parser("info", parents=[parent_parser, ksweep_parser])
    # # commands.append('info')
    
    # info_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="DELTA TREE", required=True, help="filepath to a pickle produced by the tree command. Tree nodes will be updated to hold additional sketches as needed to perform info commands selected.")

    # info_parser.add_argument("-s", "--tag", dest="tag", help="tagname used to label outputfiles, default to original tag used to create input tree",  metavar="PREFIX TAG", type=str, required=False)

    # info_parser.add_argument("-o", "--outdir", dest="outdir", default=os.getcwd(), type=str, help="directory to write the output tables.", metavar="OUTPUT DIR")

    # info_parser.add_argument("-l", "--label", dest="label", default="", help="NOT IMPLEMENTED Label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list).", required=False, metavar="SUFFIX TAG")

    # info_parser.set_defaults(func=info_command)

   # Make parser for "dand_cmd.py kij ..."
    kij_parser = subparsers.add_parser("kij", help="K Independent Jaccard. If a subset of fastas is not provided, matrix will include all inputs used to generate the delta tree using the `tree` command. NOTE: Options used during creation of delta tree will be used (e.g. exact/estimate, genome directory, species tag name.)", parents=[parent_parser, ksweep_parser])
    
    commands.append('kij')

    kij_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="DELTA TREE", required=True, help="filepath to a pickle produced by the tree command")
    kij_parser.add_argument("-s", "--tag", dest="tag", help="tagname used to label outputfiles, default to original tag used to create input tree",  metavar="PREFIX TAG", type=str, required=False)

    kij_parser.add_argument("-f", "--fastas", dest="flist_loc", default=None, type=str, metavar="FILEPATH", help="filepath to a subset of fasta files from the original tree which should be analyzed.")

    kij_parser.add_argument("-o", "--outdir", dest="outdir", default=os.getcwd(), type=str, help="directory to write the output tables.", metavar="OUTPUT DIR")

    kij_parser.add_argument("-l", "--label", dest="label", metavar="SUFFIX TAG", default="", help="NOT IMPLEMENTED Label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list).", required=False)

    kij_parser.add_argument("--afproject", dest="afproject", default=False, action="store_true", help="Indicate whether intermediate pickle of tuples should be created for use with helpers/afproject.py.")

    kij_parser.add_argument("--jaccard", dest="jaccard", default=False, action="store_true", help="Indicate whether to include output for standard jaccard difference for indicated ks")

    kij_parser.set_defaults(func=kij_command)

    return parser, commands


# def parse_args():
#     parser = argparse.ArgumentParser(
#         description="DandD: Efficient measurement of sequence growth and similarity"
#     )
    
#     # Global options
#     parser.add_argument('--debug', action='store_true',
#                        help='Print commands for each call to Dashing or KMC')
#     parser.add_argument('--fast', action='store_true',
#                        help='Avoid writing intermediate reference files')
#     parser.add_argument('--safe', action='store_true',
#                        help='Double check all sketch hashes')
#     parser.add_argument('--verbose', action='store_true',
#                        help='Show progress messages')

#     subparsers = parser.add_subparsers(dest='command')

#     # Tree command
#     tree_parser = subparsers.add_parser('tree')
#     tree_parser.add_argument('--fastas', '-f', help='File containing paths to fastas')
#     tree_parser.add_argument('--datadir', help='Directory containing fastas')
#     tree_parser.add_argument('--outdir', '-o', default='.',
#                             help='Output directory (default: current)')
#     # Add other tree arguments...

#     # Progressive command 
#     prog_parser = subparsers.add_parser('progressive')
#     prog_parser.add_argument('--dtree', '-d', required=True,
#                             help='Delta tree pickle from tree command')
#     # Add other progressive arguments...

#     # KIJ command
#     kij_parser = subparsers.add_parser('kij')
#     kij_parser.add_argument('--dtree', '-d', required=True,
#                            help='Delta tree pickle from tree command')
#     # Add other kij arguments...

#     return parser.parse_args()

def main_recommended():
    args = parse_args()
    
    if not args.command:
        print("Error: Please specify a command (tree, progressive, or kij)")
        sys.exit(1)
        
    dandd = DandD(
        debug=args.debug,
        fast=args.fast,
        safe=args.safe,
        verbose=args.verbose
    )
    
    if args.command == 'tree':
        dandd.run_tree(args)
    elif args.command == 'progressive':
        dandd.run_progressive(args)
    elif args.command == 'kij':
        dandd.run_kij(args)

def main():
    parser, commands = parse_arguments()
    args = parser.parse_args(sys.argv[1:])
    # if len(sys.argv) < 2:
    #     print('Must specify a command: ' + str(commands), file=sys.stderr)
    #     return 1
    args.func(args)
    

if __name__ == '__main__':
    main() 