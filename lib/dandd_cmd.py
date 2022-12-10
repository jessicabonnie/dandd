import argparse
import huffman_dandd
import sys
import pickle
import os
from tabulate import tabulate

# subcommand help: https://stackoverflow.com/questions/362426/implementing-a-command-action-parameter-style-command-line-interfaces

def tree_command(args):
    if not args.sketchdir:
        args.sketchdir=os.path.join(args.outdir,args.tag,"sketchdb")
        os.makedirs(args.sketchdir, exist_ok=True)
    tool="dashing"
    if args.exact:
        tool='kmc'
        
    os.makedirs(args.outdir, exist_ok=True)
    dtree = huffman_dandd.create_delta_tree(tag=args.tag, genomedir=args.genomedir, sketchdir=args.sketchdir, kstart=args.kstart, nchildren=args.nchildren, registers=args.registers, flist_loc=args.flist_loc, canonicalize=args.canonicalize, tool=tool)

    dtree.save(outdir=args.outdir, tag=args.tag, label=args.label)

def progressive_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
    results=dtree.progressive_wrapper(flist_loc=args.flist_loc, count=args.norderings, ordering_file=args.ordering_file, step=args.step)
    if args.outfile:
        results.to_csv(args.outfile)
    else:
        print(tabulate(results, headers=list(results.columns)))

def abba_command(args):
    dtree = pickle.load(open(args.delta_tree, "rb"))
    # results=dtree.progressive_wrapper(flist_loc=args.flist_loc, count=args.norderings, ordering_file=args.ordering_file, step=args.step)
    # if args.outfile:
    #     results.to_csv(args.outfile)
    # else:
    #     print(tabulate(results, headers=list(results.columns)))
    

def parse_arguments():
    parser = argparse.ArgumentParser(description='program to explore delta values for a set of fasta files')
    
    # Arguments for top-level

    parser.add_argument("-v", "--verbose", action="count", default=0)
    

    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands',help='additional help')

    # Make parser for "dand_cmd.py tree ..."
    tree_parser = subparsers.add_parser("tree")

    tree_parser.add_argument("-s", "--species", dest="tag", help="species tagname used to locate data directory",  metavar="SPECIESNAME", type=str, required=True)
    # choices=['ecoli', 'salmonella', 'human', 'HVSVC2','HVSVC2_snv', 'HVSVC2_snv_indel','HVSVC2_snv_sv', 'bds']
    tree_parser.add_argument("-x", "--exact", dest="exact", help="instead of estimating, count kmers using kmc3", default=False, action="store_true", required=False)

    tree_parser.add_argument("-g", "--genomedir", dest="genomedir", default='/scratch16/blangme2/jessica/data', help="data directory containing the fasta files -- all will be included if --fastas is not used", type=str, metavar="FASTADIR")

    tree_parser.add_argument("-o", "--out", dest="outdir", default=os.getcwd(), help="top level output directory that will contain the species directory after running", type=str, metavar="OUTDIRPATH")

    tree_parser.add_argument("-c", "--sketchdir", dest="sketchdir", default=None, help="sketch directory for species", type=str, metavar="SKETCHDIR")

    tree_parser.add_argument("-k", "--kstart", dest="kstart", default=12, help="kmer length at which to start the search for delta (different species have different optimal k values)", type=int, metavar="INT")
    
    tree_parser.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout. NOT IMPLEMENTED.")

    tree_parser.add_argument("-f", "--fastas", dest="flist_loc", metavar="FILE", type=str, default=None, help="filepath to a subset of fasta files to use in the species directory -- no title, one per line")

    tree_parser.add_argument("-l", "--label", dest="label", metavar="STRING", default=None, help="label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list). NOT IMPLEMENTED", required=False)

    tree_parser.add_argument("-n", "--nchildren", dest="nchildren", metavar="INTEGER", type=int, default=None, help="number of children for each node in the delta tree -- default is to create a tree of only 2 levels with all individual sketches as the children of the root node.")

    tree_parser.add_argument("-r", "--registers", dest="registers", metavar="INTEGER", default=20, help="number of registers to use during sketching")

    tree_parser.add_argument("-C", "--no-canon", action="store_false", default=True,  dest="canonicalize", help="instruct dashing to use non-canonicalized kmers")

    tree_parser.set_defaults(func=tree_command)


    # Make parser for "dand_cmd.py progressive ..."
    progressive_parser = subparsers.add_parser("progressive", help="Measure Delta as each individual fasta is added to the set. If a specific ordering is not provided, a set of random orderings can be generated. NOTE: Options used during creation of delta tree will be used (e.g. exact/estimate, genome directory, species tag name.")
    
    progressive_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="TREEPICKLE", required=True, help="filepath to a pickle produced by the tree command")

    progressive_parser.add_argument("-r", "--orderings", dest="ordering_file", metavar="ORDERING PICKLE", type=str, default=None, help="filepath to a pickle of orderings if different from default named using tag")

    progressive_parser.add_argument("-f", "--fastas", dest="flist_loc", default=None, type=str, metavar="FASTA LIST FILE", help="filepath to a subset of fasta files from the original tree which should be analyzed using progressive union. When count is not provided, the ordering in the file will be used for a single progression. The ordering will not be added to the ordering pickle.")

    progressive_parser.add_argument("-n", "--norderings", dest="norderings", default=None, type=int, help="number of random orderings to explore. If not provided, the orderings stored in the ordering pickle will be used. If that file does not exist / is not provided, program will terminate.")

    progressive_parser.add_argument("-o", "--outfile", dest="outfile", default=None, type=str, help="path to write the output table. If path not provided, table will be printed to standard out.")

    progressive_parser.add_argument("--step", dest="step", default=1, type=int, help="Number of sketches to include in each progression. Mostly used for a single ordered progression.")

    # progressive_parser.add_argument("-e", "--exhaustive", action="store_true", dest="exhaustive", default=False, help="instead of a randome ordering, all possible permutations of the indicated fastas will be run. Warning, this should be restricted to a small subset unless space and time considerations are managed.")

    progressive_parser.set_defaults(func=progressive_command)

    # Make parser for "dand_cmd.py info ..."
    info_parser = subparsers.add_parser("info")
    
    # Make parser for "dand_cmd.py exhaustive ..."
    abba_parser = subparsers.add_parser("abba", help='NOT IMPLEMENTED. "A before B, B before A" runs runs all permutations of orderings of subsets where both fasta A and fasta B are present -- with analysis comparing those where (1) fasta A preceeds B in the ordering; (2) fasta B preceeds A in the ordering. NOTE: A and B should both be present in the provided tree.')

    abba_parser.add_argument("-d", "--dtree", dest="delta_tree", metavar="TREEPICKLE", required=True, help="filepath to a pickle produced by the tree command")

    abba_parser.add_argument("-A", "--fastaA", dest="fastaA", metavar="FASTA_A_PATH", type=str, default=None, help="filepath of fasta A")

    abba_parser.add_argument("-B", "--fastaB", dest="fastaB", metavar="FASTA_B_PATH", type=str, default=None, help="filepath of fasta B")

    abba_parser.set_defaults(func=abba_command)

    return parser



def main():
    parser=parse_arguments()
    args = parser.parse_args(sys.argv[1:])
    args.func(args)



if __name__ == '__main__':
    main()