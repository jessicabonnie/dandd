from argparse import ArgumentParser
import huffman_dandd

# subcommand help: https://stackoverflow.com/questions/362426/implementing-a-command-action-parameter-style-command-line-interfaces

def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument("-s", "--species", dest="tag", help="species tagname used to locate data directory", choices=['ecoli', 'salmonella', 'human', 'HVSVC2','HVSVC2_snv', 'HVSVC2_snv_indel','HVSVC2_snv_sv'], metavar="SPECIESNAME", type=str, required=True)

    parser.add_argument("-g", "--genomedir", dest="genomedir", default='/home/jbonnie1/scr16_blangme2/jessica/data',help="data directory containing the species subdirectories of fasta files", type=str, metavar="DATADIR")

    parser.add_argument("-c", "--sketchdir", dest="sketchdir", default='/home/jbonnie1/scr16_blangme2/jessica/dandd/progressive_union/sketches',help="parent directory containing the species subdirectories of sketches", type=str, metavar="DATADIR")

    parser.add_argument("-k", "--kstart", dest="kstart", default=12, help="kmer length at which to start the search for delta (different species have different optimal k values)", type=int, metavar="INT")
    
    parser.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")

    parser.add_argument("-f", "--fastas", dest="flist_loc", metavar="FILE", default=None, help="filepath to a subset of fasta files to use in the species directory -- no title, one per line")

    parser.add_argument("-l", "--label", dest="label", metavar="STRING", default=None, help="label to use in result file names -- to distinguish it from others (e.g. to indicate a particular input file list)")

    args = parser.parse_args()
    return args

#tag: str, genomedir: str, sketchdir: str, kstart: int

def main():
    args = parse_arguments()
    dtree = huffman_dandd.create_delta_tree(tag=args.tag, genomedir=args.genomedir, sketchdir=args.sketchdir, kstart=args.kstart, flist_loc=args.flist_loc)
    huffman_dandd.save_dtree(dtree=dtree, outloc=args.sketchdir, tag=args.tag, label=args.label)



if __name__ == '__main__':
    main()