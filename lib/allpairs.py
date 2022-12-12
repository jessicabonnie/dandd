#!/usr/bin/env python3

"""
Driver for experiments that perform all-pairs-k-independent-jaccard (KIJ)
comparisons across a set of FASTA files.

Example invocation:

./allpairs.py --ref-tree assembled-fish_mito-ref.newick

TODO: allow one method for finding ks and another for sweeping for KIJs
"""

import os
import math
import shutil
import json
import sys
import subprocess
import itertools
from subprocess import PIPE
import argparse
import multiprocessing as mp
# you may have to pip or conda install this
from ete3 import Tree


# assume they're in PATH
bin_dict = {'dashing': 'dashing', 'dashing2': 'dashing2', 'kmc': 'kmc'}


def dashing_card_cmd(k, nestimators, estimator_width, extra, protein, inputs):
    sketch_cmd = '%s sketch -k %d -S %d %s %s %s' % (bin_dict['dashing'], k, int(math.log2(nestimators)), extra, protein, inputs)
    # dashing card -k 10 -S 18 --min-count 2 --cache-sketches camaldule.fasta
    return '%s hll -k %d -S %d %s %s %s' % (bin_dict['dashing'], k, int(math.log2(nestimators)), extra, protein, inputs)


def kmc3_card_cmd(k, nestimators, estimator_width, extra, protein, inp):
    output_pre = inp + '.kmc_pre'
    test_str = "test -f %s" % output_pre
    kmc_str = "kmc -v -k%d -fm -ci1 -cs2 %s %s /tmp/" % (k, inp, inp)
    kmctools_str = "kmc_tools info %s | head -n 2 | tail -n 1 | awk '{print $NF}'"
    cmd = '(%s || %s) && %s' % (test_str, kmc_str, kmctools_str)
    return cmd


def dashing_card(cmd):
    """
    Run "dashing hll" and return the cardinality that it computes.
    """
    p = subprocess.Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    out, err = out.decode(), err.decode()
    if p.returncode != 0:
        raise RuntimeError('Command "%s" returned %d, out="%s", err="%s"' % (cmd, p.returncode, out, err))
    out = out.rstrip()
    card_result = out.split()[-1]
    return float(card_result)


def kmc3_card(cmd):
    """
    Run "dashing hll" and return the cardinality that it computes.
    """
    p = subprocess.Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    out, err = out.decode(), err.decode()
    if p.returncode != 0:
        raise RuntimeError('Command "%s" returned %d, out="%s", err="%s"' % (cmd, p.returncode, out, err))
    return float(out.rstrip())


def dashing2_card_cmd(_1, _2, _3, _4):
    raise RuntimeError('No card_cmd function yet for dashing2')


def dashing2_card(_cmd):
    raise RuntimeError('No card function yet for dashing2')


card_command_dict = {'dashing': dashing_card_cmd,
                     'dashing2': dashing2_card_cmd,
                     'kmc': kmc3_card_cmd}

card_dict = {'dashing': dashing_card,
             'dashing2': dashing2_card,
             'kmc': kmc3_card}


def card_cmd(tool, k, nestimators, estimator_width, extra, protein, inputs):
    return card_command_dict[tool](k, nestimators, estimator_width, extra, protein, inputs)


def mp_runner(tup):
    """
    Wrapper function that runs the tool & parameters indicated by tuple
    argument.  For use with multiprocessing.Pool.
    """
    tool, name1, name2, k, cmd = tup
    return tool, name1, name2, k, card_dict[tool](cmd)


def delta_summarize(results):
    """
    Given cardinalities for the various ks, construct the deltas
    """
    best_dkk = {}
    for tool, name1, name2, k, card in results:
        key = (tool, name1, name2)
        dk_over_k = card/k
        if key not in best_dkk or dk_over_k > best_dkk[key][0]:
            best_dkk[key] = (dk_over_k, card, k)
    delta_summary = []
    for k, v in sorted(best_dkk.items()):
        tool, name1, name2 = k
        dk_over_k, card, k = v
        delta_summary.append((tool, name1, name2, dk_over_k, card, k))
    return delta_summary


def kij_summarize(delta_summary):
    """
    Given deltas for the marginals and the pairs, construct the all pairwise
    KIJs.  Returned is a list of tuples, with each tuple of this form:

    (tool, name1, name2, k, j, k1, k2, k12)

    Since kij_summarize inherently uses 3 separate argmax k's (2 for marginals,
    1 for union), the 'k1', 'k2' and 'k12' fields are used to hold these ks and
    the 'k' field is set to 0, indicating we're not using a single fixed k.
    """
    marginals, pairs = {}, {}
    tool = None
    for _tool, name1, name2, dk_over_k, card, k in delta_summary:
        tool = _tool
        if name1 == name2:
            marginals[name1] = (dk_over_k, k)
        else:
            pairs[(name1, name2)] = pairs[(name2, name1)] = (dk_over_k, k)
    keys = marginals.keys()
    kij_summary = []
    for key1, key2 in itertools.combinations(keys, 2):
        dkk1, k1 = marginals[key1]
        dkk2, k2 = marginals[key2]
        dkk12, k12 = pairs[(key1, key2)]
        kij = (dkk1 + dkk2 - dkk12)/dkk12
        kij_summary.append((tool, key1, key2, 0, kij, k1, k2, k12))
    return kij_summary


def j_summarize(results, target_k):
    """
    Given fixed-k cardinalities for the marginals and the pairs, construct the
    all pairwise Js.  Returned is a list of tuples, with each tuple of this
    form:

    (tool, name1, name2, k, j, k1, k2, k12)

    Since j_summarize inherently uses a single fixed k, the 'k' field will be
    set to this positive value of k, and k1, k2 and k12 will be None.
    """
    marginals, pairs = {}, {}
    tool = None
    for _tool, name1, name2, k, card in results:
        tool = _tool
        if k != target_k:
            continue
        if name1 == name2:
            marginals[name1] = card
        else:
            pairs[(name1, name2)] = pairs[(name2, name1)] = card
    keys = marginals.keys()
    j_summary = []
    for key1, key2 in itertools.combinations(keys, 2):
        card1, card2 = marginals[key1], marginals[key2]
        card12 = pairs[(key1, key2)]
        j = (card1 + card2 - card12)/card12
        j_summary.append((tool, key1, key2, target_k, j, None, None, None))
    return j_summary


def summ_to_phylip(summ, seqid_to_treid, phylip_fn, convert_to_ani=False):
    recs = {}
    names = set()
    assert len(summ) > 0
    for _, name1, name2, k, j, k1, k2, k12 in summ:
        names.add(name1)
        names.add(name2)
        if convert_to_ani:
            k_for_mash_distance = k
            if k == 0:
                assert k1 is not None
                assert k2 is not None
                assert k12 is not None
                k_for_mash_distance = max(k1, k2, k12)
            recs[(name1, name2)] = recs[(name2, name1)] = mash_distance(j, k_for_mash_distance)
        else:
            recs[(name1, name2)] = recs[(name2, name1)] = 1-j
    names = sorted(list(names))
    assert len(names) == len(seqid_to_treid), (len(names), len(seqid_to_treid))
    with open(phylip_fn, 'wt') as fh:
        print(str(len(names)), file=fh)
        for i, name1 in enumerate(names):
            row = [name1]
            for name2 in names[:i]:
                row.append(recs[(name1, name2)])
            print(' '.join(map(str, row)), file=fh)


def mash_distance(j, k):
    """
    This is a distance already; no need to invert the sense later
    """
    try:
        j = max(j, sys.float_info.epsilon)
        return -(math.log(2.*j/(1.+j)))/k
    except ValueError:
        raise RuntimeError('Could not compute mash distance for k=%d, j=%f' % (k, j))


def run_fneighbor(phylip_fn, tree_fn_base, seqid_to_treid, ref_tree_fn,
                  rename_seqs=False, fneighbor_exe=None, nreplicates=1):
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
        seed = i * 2 - 1
        fneighbor_fn = '%s.%d.fneighbor' % (tree_fn_base, seed)
        tree_fn = '%s.%d.newick' % (tree_fn_base, seed)
        cmd = '%s -datafile %s -outfile %s -outtreefile %s -matrixtype l -jumble -seed %d' % \
              (fneighbor_exe, phylip_fn, fneighbor_fn, tree_fn + '.tmp', seed)
        ret = os.system(cmd)
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


def rename_seqids_in_tree(orig_tree, seqid_to_treid):
    output_tree = ''
    i = 0
    while i < len(orig_tree):
        done = False
        for k in seqid_to_treid.keys():
            if orig_tree[i:].startswith(k):
                output_tree += seqid_to_treid[k]
                i += len(k)
                done = True
                break
        if done:
            continue
        output_tree += orig_tree[i]
        i += 1
    return output_tree


_default_klist = ','.join(map(str, range(2, 100)))


def go():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", type=str, default='kmc', help="which tool to use")
    parser.add_argument("--dataset", type=str, default='dataset.json', help="AFproject dataset json file")
    parser.add_argument("--write-commands", type=str, default='commands.txt', help="write commands here")
    parser.add_argument("--card-results", type=str, default='card.tsv', help="write raw cardinalities results here")
    parser.add_argument("--delta-results", type=str, default='delta.tsv', help="write delta results here")
    parser.add_argument("--j-results", type=str, default='sim.tsv', help="write final all-pairs 1-minus-Js and "
                                                                         "1-minus-KIJs to filenames like this")
    parser.add_argument("--j-results-phylip", type=str, default='sim.phylip', help="write final all-pairs "
                                                                                   "1-minus-Js and 1-minus-KIJs "
                                                                                   "here, PHYLIP format")
    parser.add_argument("--j-tree", type=str, default='kij.newick', help="write fneighbor-built 1-minus-Js and "
                                                                         "1-minus-KIJ trees to filenames like this")
    parser.add_argument("--j-tree-compare", type=str, default='kij.stats', help="write 1-minus-Js and 1-kij tree "
                                                                                "stats to filenames like this")
    parser.add_argument("--ani-results", type=str, default='ani.tsv', help="write final all-pairs 1-minus-ANIs here")
    parser.add_argument("--ani-results-phylip", type=str, default='ani.phylip', help="write all-pairs 1-minus-ANIs "
                                                                                     "here, PHYLIP format")
    parser.add_argument("--ani-tree", type=str, default='ani.newick', help="write fneighbor-built 1-minus-ANI tree here")
    parser.add_argument("--ani-tree-compare", type=str, default='ani.stats', help="write 1-ani tree comparison stats")
    parser.add_argument("--ref-tree", type=str, default='tree.newick', help="reference tree, newick format, for input")
    parser.add_argument("--fneighbor", type=str, default=None, help="path to fneighbor executable; "
                                                                    "default: look in PATH")
    parser.add_argument("--fneighbor-replicates", type=int, default=1, help="# times to invoke fneighbor with new "
                                                                             "random seeds")
    parser.add_argument("--replicates", type=int, default=1, help="# times to invoke the whole shebang")
    parser.add_argument("--rename-seqs", type=bool, default=False, help="rename to long names before tree compare")
    parser.add_argument("--extra", type=str, default='', help="extra arguments for sketching tool")
    parser.add_argument("--nest", type=int, default=262144, help="# estimators (power of 2 required for dashing)")
    parser.add_argument("--bitsper", type=int, default=8, help="bits per estimator (default 8; 8 required for dashing)")
    parser.add_argument("--klist", type=str, default=_default_klist, help="ks to try")
    parser.add_argument("--cpu", type=int, default=-1,
                        help="number of simultaneous cores to use for sketch/cardinality commands")
    args = parser.parse_args()
    print('Performing run with name "%s"' % args.name)
    if os.path.exists(args.name):
        raise RuntimeError('Output directory with name "%s" already exists' % args.name)
    os.makedirs(args.name)
    if not os.path.exists(args.dataset):
        raise RuntimeError('No dataset file "%s"' % args.dataset)
    with open(args.dataset, 'rt') as fh:
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
        inp += '.fasta'
        if not os.path.exists(inp):
            raise RuntimeError('Input path does not exist: "%s"' % inp)
        inputs.append(inp)
        inp_base = os.path.basename(inp)
        assert inp_base.count('.') == 1
        input_names.append(inp_base.split('.')[0])
    # Compose a bunch of dashing hll commands
    commands = []
    for k in map(int, args.klist.split(',')):
        for i in range(len(inputs)):
            i1, i1_name = inputs[i], input_names[i]
            commands.append((args.tool, i1_name, i1_name, k,
                             card_cmd(args.tool, k, args.nest, args.bitsper, args.extra, '', i1)))
            for j in range(i+1, len(inputs)):
                i2, i2_name = inputs[j], input_names[j]
                commands.append((args.tool, i1_name, i2_name, k,
                                 card_cmd(args.tool, k, args.nest, args.bitsper, args.extra, '', i1 + ' ' + i2)))
    # write all the commands to a file
    write_commands_fn = os.path.join(args.name, args.write_commands)
    with open(write_commands_fn, 'wt') as fh:
        print('\n'.join(map(str, commands)), file=fh)
    # Run the commands using Python multithreading
    cpu = args.cpu
    if cpu < 0:
        cpu = mp.cpu_count()
    with mp.Pool(cpu) as p:
        results = p.map(mp_runner, commands)
    with open(args.card_results, 'wt') as fh:
        print('\t'.join(['tool', 'name1', 'name2', 'k', 'card']), file=fh)
        for tool, name1, name2, k, card in results:
            print('\t'.join([tool, name1, name2, str(k), str(card)]), file=fh)
    # Pick out maximal dk/ks
    dsumm = delta_summarize(results)
    with open(args.delta_results, 'wt') as fh:
        print('\t'.join(['tool', 'name1', 'name2', 'delta', 'card', 'k']), file=fh)
        for tool, name1, name2, delta, card, k in dsumm:
            print('\t'.join([tool, name1, name2, str(k), str(card), str(delta)]), file=fh)
    # Final combination that combines deltas for marginals and pairs
    j_and_kij_summ = kij_summarize(dsumm)
    # Records are of the form (tool, name1, name2, k, j, k1, k2, k12)
    for k in map(int, args.klist.split(',')):
        j_and_kij_summ += j_summarize(results, k)
    for k in [0] + list(map(int, args.klist.split(','))):
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
            stats = run_fneighbor(phylip_fn, j_tree_fn, seqid_to_treid, args.ref_tree,
                                  rename_seqs=args.rename_seqs,
                                  fneighbor_exe=args.fneighbor,
                                  nreplicates=args.fneighbor_replicates)
            with open(j_stats_fn, 'wt') as fh:
                for i, (nrf, rf, max_rf) in enumerate(stats):
                    print('\t'.join(map(str, ['j', k, k1, k2, k12, nrf, rf, max_rf, i])), file=fh)

        # write Phylip output for ANI
        if args.ani_results_phylip:
            fn = _customize_fn(args.ani_results_phylip)
            summ_to_phylip(summ, seqid_to_treid, fn, convert_to_ani=True)
            assert os.path.exists(fn)

        # write fneighbor tree output for ANI
        ani_tree_fn = _customize_fn(args.ani_tree)
        if args.ani_tree:
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


if __name__ == '__main__':
    go()
