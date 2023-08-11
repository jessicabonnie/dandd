import os, sys, pickle
sys.path.insert(0, '../lib')
import shutil
import json
from allpairs import summ_to_phylip, rename_seqids_in_tree
# you may have to pip or conda install this
from ete3 import Tree
import argparse




def run_fneighbor(phylip_fn, tree_fn_base, seqid_to_treid, ref_tree_fn,
                  rename_seqs=False, fneighbor_exe=None, nreplicates=1, debug=False):
    """
    mamba install embassy-phylip
    fneighbor -datafile input.phylip -outfile output.fneighbor \
        -outtreefile output.newick -matrixtype l -jumble -seed N
    """
    if fneighbor_exe is None:
        fneighbor_exe = shutil.which("fneighbor")
    if fneighbor_exe is None:
        raise RuntimeError('Could not find fneighbor binary in path')
    results = []
    for i in range(nreplicates):
        seed = i * 2 + 1
        fneighbor_fn = '%s.%d.fneighbor' % (tree_fn_base, seed)
        tree_fn = '%s.%d.newick' % (tree_fn_base, seed)
        cmd = '%s -datafile %s -outfile %s -outtreefile %s -matrixtype l -jumble -seed %d' % \
              (fneighbor_exe, phylip_fn, fneighbor_fn, tree_fn + '.tmp', seed)
        ret = os.system(cmd)
        if debug:
            print(cmd)
        if ret != 0:
            raise RuntimeError('Non-zero return code (%d) from command "%s"' % (ret, cmd))
        assert os.path.exists(fneighbor_fn)
        assert os.path.exists(tree_fn + '.tmp')
        if rename_seqs:
            with open(tree_fn, 'wt') as tree_fh:
                orig_tree = open(tree_fn + '.tmp', 'rt').read()
                output_tree = rename_seqids_in_tree(orig_tree, seqid_to_treid)
                print(output_tree, file=tree_fh)
        else:
            shutil.copyfile(tree_fn + '.tmp', tree_fn)
        nrf, rf, max_rf = ete_compare(tree_fn, ref_tree_fn)
        results.append([nrf, rf, max_rf])
    return results


def ete_compare(usr_tree_str, ref_tree_str, outgroup_id=None):
    """
    This code was adapted from https://github.com/afproject-org/afproject/blob/master/libs/treeutils.py
    Available under a Mozilla Public License
    The code is used in publicly available webservice AFproject (http://afproject.org).
    """
    qt = Tree(open(usr_tree_str, 'rt').read())
    rt = Tree(open(ref_tree_str, 'rt').read())
    if outgroup_id:
        qt.set_outgroup(outgroup_id)
        rt.set_outgroup(outgroup_id)
    res = qt.compare(rt, unrooted=False if outgroup_id else True)
    rf = res["rf"]
    max_rf = res["max_rf"]
    nrf = res["norm_rf"]
    assert len(res["common_edges"]) > 0
    return nrf, rf, max_rf


def phylip_command(args):

    if args.phylo_ref and args.ref_tree:
        j_and_kij_summ = pickle.load(open(args.aftuples, "rb"))
        upperout=os.path.dirname(args.outfile)
        outbase=os.path.basename(args.outfile)
        dir_dict=dict([(x, os.path.join(upperout,x)) for x in ["phylip","newick", "stats", "fneighbor"]])
        #dir_dict.update({ x : os.path.join(upperout,x) } for x in ["phylip","newick", "stats", "fneighbor"])
        [os.makedirs(x,exist_ok=True) for x in dir_dict.values()]
        args.j_results_phylip=os.path.join(dir_dict["phylip"], outbase + '.sim.phylip')
        args.j_tree=os.path.join(dir_dict["newick"], outbase + '.newick')
        #args.j_results=args.outfile + '.sim.tsv'
        args.j_tree_compare=os.path.join(dir_dict["stats"], outbase + '.stats')
        args.replicates=1
        args.fneighbor_replicates=1
        args.nest=262144
        args.bitsper=8
        args.cpu=-1
        args.fneighbor='fneighbor'
        args.datadir=os.path.dirname(args.dtree.fastas[0])
        args.ani_results_phylip=None
        args.ani_tree=None

        phylogeny_comparison(args, j_and_kij_summ=j_and_kij_summ)

        


def phylogeny_comparison(args, j_and_kij_summ):
    '''Create necessary files for comparison using AFproject.
    NOTE: very brittle in terms of expectation of relationship of file names to sample
    '''
    if not os.path.exists(args.phylo_ref):
        raise RuntimeError('No dataset file "%s"' % args.dataset)
    with open(args.phylo_ref, 'rt') as fh:
        data = json.load(fh)
    seqid_to_treid, treid_to_seqid = {}, {}
    assert len(data['treids']) == 0 or len(data['treids']) == len(data['seqids'])
    if len(data['treids']) == 0:
        for sid in data['seqids']:
            seqid_to_treid[sid] = sid
    else:
        for sid, tid in zip(data['seqids'], data['treids']):
            seqid_to_treid[sid] = tid
            treid_to_seqid[tid] = sid
    assert len(seqid_to_treid) > 0
    inputs, input_names = [], []
    for i, inp in enumerate(data['seqids']):
        inp = os.path.join(args.datadir,inp+'.fasta')
        if not os.path.exists(inp):
            raise RuntimeError('Input path does not exist: "%s"' % inp)
        inputs.append(inp)
        inp_base = os.path.basename(inp)
        assert inp_base.count('.') == 1
        input_names.append(inp_base.split('.')[0])
    for k in [0] + list(map(int, range(int(args.mink), int(args.maxk)))):
        summ = list(filter(lambda x: x[3] == k, j_and_kij_summ))
        k1, k2, k12 = summ[0][-3:]
        assert len(summ) > 0
        name = 'kij' if k == 0 else 'k%d' % k

        def _customize_fn(_fn):
            toks = _fn.split('.')
            return '.'.join(toks[:-1] + [name] + [toks[-1]])

        # write Phylip output for KIJ/J
        if args.j_results_phylip:
            fn = _customize_fn(args.j_results_phylip)
            summ_to_phylip(summ, seqid_to_treid, fn)
            assert os.path.exists(fn)

        # write fneighbor tree output for 1-Jaccard / 1-KIJ
        j_tree_fn = _customize_fn(args.j_tree)
        if args.j_tree:
            assert args.j_tree_compare is not None
            j_stats_fn = _customize_fn(args.j_tree_compare)
            phylip_fn = _customize_fn(args.j_results_phylip)
            stats = run_fneighbor(phylip_fn, j_tree_fn, seqid_to_treid, args.ref_tree, fneighbor_exe=args.fneighbor, nreplicates=args.fneighbor_replicates)
            with open(j_stats_fn, 'wt') as fh:
                for i, (nrf, rf, max_rf) in enumerate(stats):
                    print('\t'.join(map(str, ['j', k, k1, k2, k12, nrf, rf, max_rf, i])), file=fh)
        
         # write Phylip output for ANI
        if args.ani_results_phylip:
            fn = _customize_fn(args.ani_results_phylip)
            summ_to_phylip(summ, seqid_to_treid, fn, convert_to_ani=True)
            assert os.path.exists(fn)

        # write fneighbor tree output for ANI
        if args.ani_tree:
            ani_tree_fn = _customize_fn(args.ani_tree)
            assert args.ani_tree_compare is not None
            ani_stats_fn = _customize_fn(args.ani_tree_compare)
            phylip_fn = _customize_fn(args.ani_results_phylip)
            stats = run_fneighbor(phylip_fn, ani_tree_fn, seqid_to_treid, args.ref_tree,
                                  rename_seqs=args.rename_seqs,
                                  fneighbor_exe=args.fneighbor,
                                  nreplicates=args.fneighbor_replicates)
            with open(ani_stats_fn, 'wt') as fh:
                for i, (nrf, rf, max_rf) in enumerate(stats):
                    print('\t'.join(map(str, ['ani', k, k1, k2, k12, nrf, rf, max_rf, i])), file=fh)




def main():
    parser = argparse.ArgumentParser(description='program to compare output from dandd kij with AFproject benchmarks')
    
    parser.add_argument("-o", "--outdir", help="Output directory to create subdirectories. Default CWD.", default=os.getcwd(), dest="outdir")
    
    parser.add_argument("-f", "--aftuples", required=True, dest="aftuples", help="Pickle output from dandd kij --afproject command.")
    
    parser.add_argument("-d", "--dtree_loc", required=True, dest="dtree_loc", help="Pickled Delta Tree used during prior dandd kij command.")
    
    parser.add_argument("-t","--phylo-tree", dest="ref_tree", default=None, type=str, help="path to the phylogenic reference tree from AF project. Must match --phylo-ref.", required=True)

    parser.add_argument("-r","--phylo-ref", dest="phylo_ref", default=None, type=str, help="path to the reference dataset from AF project. Must match --phylo-tree.", required=True)

    parser.add_argument("--rename-seqs","--rename", dest="rename_seqs", default=False, action="store_true", help="Rename sequences to match AF project?")

    parser.add_argument("--mink", dest="mink", metavar="MINIMUM-K", required=False, help="Minimum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k", default=10)

    parser.add_argument("--maxk", dest="maxk", metavar="MAXIMUM-K", required=False, help="Maximum k to start sweep of ks for their possible deltas. Can be used to graph the argmax k", default=20)


    args = parser.parse_args(sys.argv[1:])

    args.dtree = pickle.load(open(args.dtree_loc, "rb"))
    args.outfile = args.dtree.make_prefix(tag=args.dtree.speciesinfo.tag, label=f"AF")

    phylip_command(args)




if __name__ == '__main__':
    sys.exit(main())