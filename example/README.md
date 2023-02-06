# Example: Fish Mitochondria

**Note**: this assumes that `dashing` is in your `$PATH` and that `dandd` has been made executable (`chmod +x lib/dandd`).

## Data
This example uses data from the AF project. Use the following code to download the data in this example directory:

```
cd example
wget https://afproject.org/media/genome/std/assembled/fish_mito/dataset/assembled-fish_mito.zip
unzip assembled-fish_mito.zip
```

## DandD Tree

The first step will be to create intial sketches of the fasta files and find delta for the full set. We will use "fish-mito" as the tag for filenames. It will default to writing in the current working directory. The default starting value for k is 14. This is a bit higher than the argmax-k that delta needs. DandD will discover this on its own, although we could help it by giving a starting k.

```
../lib/dandd tree --datadir ./assembled-fish_mito --tag fish-mito

```
The table of delta values for the sketches is saved to: `fish-mito_25_dashing_deltas.csv`. The `sketchdb` directory contains the dashing sketches pertaining to everything run in this output directory, as well as the `fastahex.pickle` and `sketchinfo.pickle` with metadata for DandD to track them. For a human readable table tracking the fastas and their associated sketches for each experiment see `fish-mito_25_dashing_sketchdb.txt`. (`{experiment_tag}_{#fastas}_{tool}_sketchdb.txt`)

## Children
DandD will default to unioning all of the sketches/databases at once to obtain n+1 deltas. You can also ask it to combine them in a tree-like manner by indicating the maximum number of childnodes that should be combined at each level.
```
../lib/dandd tree --datadir ./assembled-fish_mito --tag fish-mito --nchildren 2
```

## SubTree
Suppose you are only interested in a subset of the fastas for a particular experiment? You can do that too! DandD will use the sketches already produced to create this smaller tree. You will note that the only sketching/kmer counting being performed is on the union of the subset.  

```
realpath ./assembled-fish_mito/* | head -n5 > fish_subset.txt
../lib/dandd tree --fastas fish_subset.txt --tag small_fish -k 12

```
The table of delta values for the sketches is saved to: `small_fish_5_dashing_deltas.csv`. The pickle of the delta-tree object is saved to: `small_fish_5_dashing_dtree.pickle`.

## Progressive Union

We can run the progressive union on either the full set or the subset. Here we run it on the subset, since a lower number of samples means less permutations would be needed to get a look at the pattern. We start with 10 orderings. DandD will create a default file prefix from the tag name of the tree and the number of orderings, but this can also be provided.
```
../lib/dandd progressive --dtree small_fish_5_dashing_dtree.pickle --norder 10 
```

This creates `sketchdb/small_fish_5_orderings.pickle` containing the orderings. It might be the case that these orderings are desired to be used for a parallel dataset (say the same samples stratified by variant type, as in the paper). Alternatively, this allows for reproduction of these initial orderings when additional orderings are added in a later run:

```
../lib/dandd progressive --dtree small_fish_5_dashing_dtree.pickle --norder 30 
```

### K-Independent Jaccard
If we wish to calculate the K-Independent Jaccard (KIJ) for our set, we use the `kij` command. If we also want the non-k independent pairwise jaccard values we use the `--jaccard` flag and provide a window of `--mink` and `--maxk` over which we would like our non-k independent Jaccards to be calculated.

```
../lib/dandd kij --dtree small_fish_5_dashing_dtree.pickle
```
