import argparse
import huffman_dandd
import sys
import pickle
import os
# from tabulate import tabulate
import csv
from typing import List, Dict, Set, Tuple, NamedTuple

# subcommand help: https://stackoverflow.com/questions/362426/implementing-a-command-action-parameter-style-command-line-interfaces

def write_listdict_to_csv(outfile: str, listdict:List[Dict], suffix:str=""):
    writer = open(outfile+suffix, "w") if outfile is not None and outfile != '-' else sys.stdout
    fieldnames=set()
    for x in listdict:
        fieldnames.update(x.keys())
    dict_writer = csv.DictWriter(writer, fieldnames=list(fieldnames))
    dict_writer.writeheader()
    dict_writer.writerows(listdict)
    writer.close()

def insert_pre_ext(filename, string):
    toks = filename.split('.')
    return '.'.join(toks[:-1] + [string] + [toks[-1]])

def tree_command(args):
    if not args.sketchdir:
        # args.sketchdir=os.path.join(args.outdir,args.tag,"sketchdb")
        args.sketchdir=os.path.join(args.outdir,"sketchdb")
        os.makedirs(args.sketchdir, exist_ok=True)
    tool="dashing"
    if args.exact:
        tool='kmc'
        
    os.makedirs(args.outdir, exist_ok=True)
    # os.chdir(args.sketchdir)
    dtree = huffman_dandd.create_delta_tree(tag=args.tag, genomedir=args.genomedir, sketchdir=args.sketchdir, kstart=args.kstart, nchildren=args.nchildren, registers=args.registers, flist_loc=args.flist_loc, canonicalize=args.canonicalize, tool=tool, debug=args.debug, nthreads=args.nthreads)

    dtree.save(outdir=args.outdir, tag=args.tag, label=args.label)
    #subset the fastahex map to output a human readable version containing info for the sketches relevant to the tree
    # explist=[dtree.speciesinfo.sketchinfo[item] for item in list(dtree.experiment["baseset"])]
    # expmaploc=os.path.join(args.outdir,  args.tag + "_" + str(self.ngen) + "_" + self.experiment["tool"] + '_sketchdb.txt')
    # write_listdict_to_csv(outfile=args.outdir, tag=args.tag, label=args.label)

def progressive_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
        
    dtree.experiment["debug"] = args.debug
    results=dtree.progressive_wrapper(flist_loc=args.flist_loc, count=args.norderings, ordering_file=args.ordering_file, step=args.step)
    write_listdict_to_csv(outfile=args.outfile, listdict=results)
    dtree.save(outdir=os.path.split(args.outfile)[0], tag=dtree.speciesinfo.tag, label=f"progu{args.norderings}")

def abba_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
    # results=dtree.progressive_wrapper(flist_loc=args.flist_loc, count=args.norderings, ordering_file=args.ordering_file, step=args.step)
    # if args.outfile:
    #     results.to_csv(args.outfile)
    # else:
    #     print(tabulate(results, headers=list(results.columns)))

def info_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
    summary=dtree.summarize(mink=args.mink, maxk=args.maxk)
    write_listdict_to_csv(outfile=args.outfile, listdict=summary)

def kij_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
    fastas=[]
    if args.flist_loc:
        with open(args.flist_loc) as file:
            fastas = [line.strip() for line in file]

    dtree.ksweep(mink=args.mink,maxk=args.maxk)
    kij_results, j_results = dtree.pairwise_spiders(sublist=fastas, mink=args.mink, maxk=args.maxk)
    write_listdict_to_csv(outfile=args.outfile, listdict=kij_results)

    # writer = open(args.outfile, "w") if args.outfile is not None and args.outfile != '-' else sys.stdout
    # dict_writer = csv.DictWriter(writer, fieldnames=kij_results[0].keys())
    # dict_writer.writeheader()
    # dict_writer.writerows(kij_results)
    # writer.close()
    
    if args.jaccard:
        write_listdict_to_csv(outfile=args.outfile,listdict=j_results,suffix=".jaccard")
        # writer = open(args.outfile+".jaccard", "w") if args.outfile is not None and args.outfile != '-' else sys.stdout
        # dict_writer = csv.DictWriter(writer, fieldnames=j_results[0].keys())
        # dict_writer.writeheader()
        # dict_writer.writerows(j_results)
        # writer.close()

    
    

def parse_arguments():
    parser = argparse.ArgumentParser(description='program to explore delta values for a set of fasta files')
    
    # Arguments for top-level

    commands = []
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument( "--debug", action="store_true", default=False, dest="debug")

    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands',help='additional help')

    # Make parser for "dand_cmd.py tree ..."
    tree_parser = subparsers.add_parser("tree")
    commands.append('tree')
    # tree_parser.add_argument( "--debug", action="store_true", default=False, dest="debug")

    tree_parser.add_argument("-s","-t", "--tag", dest="tag", help="tagname used to label outputfiles; if datadir contains subdirectory by the same name fastas will be sourced from there",  metavar="LABELSTRING", type=str, required=False, default='dandd')
    # choices=['ecoli', 'salmonella', 'human', 'HVSVC2','HVSVC2_snv', 'HVSVC2_snv_indel','HVSVC2_snv_sv', 'bds']
    tree_parser.add_argument("-x", "--exact", dest="exact", help="instead of estimating, count kmers using kmc3", default=False, action="store_true", required=False)

    tree_parser.add_argument("-d", "--datadir", dest="genomedir", default='/scratch16/blangme2/jessica/data', help="data directory containing the fasta files -- all will be included if --fastas is not used", type=str, metavar="FASTADIR")

    tree_parser.add_argument("-o", "--out", dest="outdir", default=os.getcwd(), help="top level output directory that will contain the output files after running", type=str, metavar="OUTDIRPATH")

    tree_parser.add_argument("-c", "--sketchdir", dest="sketchdir", default=None, help="sketch directory for species", type=str, metavar="SKETCHDIR")

    tree_parser.add_argument("-k", "--kstart", dest="kstart", default=12, help="kmer length at which to start the search for delta (different species have different optimal k values)", type=int, metavar="INT")
    
    tree_parser.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout. NOT IMPLEMENTED.")

    tree_parser.add_argument("-f", "--fastas", dest="flist_loc", metavar="FILE", type=str, default=None, help="filepath to a subset of fasta files to use in the species directory -- no title, one per line")

    tree_parser.add_argument("-l", "--label", dest="label", metavar="STRING", default="", help="NOT IMPLEMENTED. label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list). NOT IMPLEMENTED", required=False)

    tree_parser.add_argument("-n", "--nchildren", dest="nchildren", metavar="INTEGER", type=int, default=None, help="number of children for each node in the delta tree -- default is to create a tree of only 2 levels with all individual sketches as the children of the root node.")

    tree_parser.add_argument("-r", "--registers", dest="registers", metavar="INTEGER", default=20, help="number of registers to use during sketching")

    tree_parser.add_argument("-e", "--nthreads", dest="nthreads", metavar="INTEGER", type=int, default=10, help="number of threads to use in calls to associated kmer counting/estimating tool")


    tree_parser.add_argument("-C", "--no-canon", action="store_false", default=True,  dest="canonicalize", help="instruct dashing to use non-canonicalized kmers")

    tree_parser.set_defaults(func=tree_command)


    # Make parser for "dand_cmd.py progressive ..."
    progressive_parser = subparsers.add_parser("progressive", help="Measure Delta as each individual fasta is added to the set. If a specific ordering is not provided, a set of random orderings can be generated. NOTE: Options used during creation of delta tree will be used (e.g. exact/estimate, genome directory, species tag name.)")
    commands.append('progressive')

    progressive_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="TREEPICKLE", required=True, help="filepath to a pickle produced by the tree command")

    progressive_parser.add_argument("-r", "--orderings", dest="ordering_file", metavar="ORDERING PICKLE", type=str, default=None, help="filepath to a pickle of orderings if different from default named using tag")

    progressive_parser.add_argument("-f", "--fastas", dest="flist_loc", default=None, type=str, metavar="FASTA LIST FILE", help="filepath to a subset of fasta files from the original tree which should be analyzed using progressive union. When count is not provided, the ordering in the file will be used for a single progression. The ordering will not be added to the ordering pickle.")

    progressive_parser.add_argument("-n", "--norderings", dest="norderings", default=0, type=int, help="number of random orderings to explore. If not provided, the orderings stored in the ordering pickle will be used. If that file does not exist / is not provided, program will terminate.")

    progressive_parser.add_argument("-o", "--outfile", dest="outfile", default=None, type=str, help="path to write the output table. If path not provided, table will be printed to standard out.")

    # progressive_parser.add_argument("-o", "--out", dest="outdir", default=os.getcwd(), help="top level output directory that will contain the output files after running", type=str, metavar="OUTDIRPATH")

    progressive_parser.add_argument("--step", dest="step", default=1, type=int, help="Number of sketches to include in each progression. Mostly used for a single ordered progression.")

    # progressive_parser.add_argument("-e", "--exhaustive", action="store_true", dest="exhaustive", default=False, help="instead of a randome ordering, all possible permutations of the indicated fastas will be run. Warning, this should be restricted to a small subset unless space and time considerations are managed.")

    progressive_parser.set_defaults(func=progressive_command)

    # Make parser for "dand_cmd.py info ..."
    info_parser = subparsers.add_parser("info")
    commands.append('info')
    
    info_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="TREEPICKLE", required=True, help="filepath to a pickle produced by the tree command. Tree nodes will be updated to hold additional sketches as needed to perform info commands selected.")
    
    info_parser.add_argument("--mink", dest="mink", metavar="MINIMUM-K", required=False, default=10, help="Minimum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    info_parser.add_argument("--maxk", dest="maxk", metavar="MAXIMUM-K", required=False, help="Maximum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    info_parser.add_argument("-o", "--outfile", dest="outfile", default=None, type=str, help="path to write the output table. If path not provided, table will be printed to standard out.")


    info_parser.set_defaults(func=info_command)
    
    # Make parser for "dand_cmd.py abba ..."
    abba_parser = subparsers.add_parser("abba", help='NOT IMPLEMENTED. "A before B, B before A" runs runs all permutations of orderings of subsets where both fasta A and fasta B are present -- with analysis comparing those where (1) fasta A preceeds B in the ordering; (2) fasta B preceeds A in the ordering. NOTE: A and B should both be present in the provided tree.')
    commands.append('abba')

    abba_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="TREEPICKLE", required=True, help="filepath to a pickle produced by the tree command")

    abba_parser.add_argument("-A", "--fastaA", dest="fastaA", metavar="FASTA_A_PATH", type=str, default=None, help="filepath of fasta A")

    abba_parser.add_argument("-B", "--fastaB", dest="fastaB", metavar="FASTA_B_PATH", type=str, default=None, help="filepath of fasta B")

    abba_parser.set_defaults(func=abba_command)

   # Make parser for "dand_cmd.py kij ..."
    kij_parser = subparsers.add_parser("kij", help="K Independent Jaccard. If a subset of fastas is not provided, matrix will include all inputs used to generate the delta tree using the `tree` command. NOTE: Options used during creation of delta tree will be used (e.g. exact/estimate, genome directory, species tag name.)")
    commands.append('kij')

    kij_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="TREEPICKLE", required=True, help="filepath to a pickle produced by the tree command")

    kij_parser.add_argument("-f", "--fastas", dest="flist_loc", default=None, type=str, metavar="FASTA LIST FILE", help="filepath to a subset of fasta files from the original tree which should be analyzed.")

    kij_parser.add_argument("-o", "--outfile", dest="outfile", default=None, type=str, help="path to write the output table. If path not provided, table will be printed to standard out.")

    # kij_parser.add_argument("-o", "--out", dest="outdir", default=os.getcwd(), help="top level output directory that will contain the output files after running", type=str, metavar="OUTDIRPATH")

    kij_parser.add_argument("--jaccard", dest="jaccard", default=False, action="store_true", help="Indicate whether to include output for standard jaccard difference for indicated ks")

    kij_parser.add_argument("--mash", dest="mash", default=False, action="store_true", help="Indicate whether to include output for mash difference for indicated ks")

    kij_parser.add_argument("--mink", dest="mink", metavar="MINIMUM-K", required=False, type=int, default=10, help="Minimum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    kij_parser.add_argument("--maxk", dest="maxk", metavar="MAXIMUM-K", required=False, type=int, default=32, help="Maximum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k")

    kij_parser.set_defaults(func=kij_command)

    return parser, commands



def main():
    parser, commands = parse_arguments()
    args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 2:
        print('Must specify a command: ' + str(commands), file=sys.stderr)
        return 1
    args.func(args)
    return 0



if __name__ == '__main__':
    sys.exit(main())