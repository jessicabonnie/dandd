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

If you would rather provide a list of global paths to the input fastas rather than a directory, the following command will produce the same results as above. If both commands are run, observe that no additional sketches need to be produced.
```
realpath ./assembled-fish_mito/* > fish_list.txt
../lib/dandd tree --fastas fish_list.txt --tag fish-mito
```

## Children
DandD will default to unioning all of the sketches/databases at once to obtain n+1 deltas. You can also ask it to combine them in a tree-like manner by indicating the maximum number of childnodes that should be combined at each level. In this example we also instruct DandD to start looking for the optimal k at `k=12`. Note that the initial sketches of the individual fastas were created during the previous example, so all that is needed is to produce the desired unions. You will observe that the output file `fish-mito_25_dashing_deltas.csv` has around twice as many lines containing delta values as in the previous command because the deltas for the additional intermediate unionsare also reported.
```
../lib/dandd tree --datadir ./assembled-fish_mito --tag fish-mito2 --nchildren 2 -k 12
```

## SubTree
Suppose you are only interested in a subset of the fastas for a particular experiment? You can do that too! DandD will use the sketches already produced to create this smaller tree. You will note that the only sketching/kmer-counting being performed is on the union of the subset, since no `--nchildren` value is indicated and all individual fasta sketches already exist at the relative values of k.  

```
realpath ./assembled-fish_mito/* | head -n5 > fish_subset.txt
../lib/dandd tree --fastas fish_subset.txt --tag small_fish -k 12

```
The table of delta values for the sketches is saved to: `small_fish_5_dashing_deltas.csv`. The pickle of the delta-tree object is saved to: `small_fish_5_dashing_dtree.pickle`. The number of input fastas is always included in these output files.

## Progressive Union

We can run the progressive union on either the full set or the subset. Here we run it on the subset, since a lower number of samples means less permutations would be needed to get a look at the pattern. We start with 10 orderings. DandD will create a default file prefix from the tag name of the tree and the number of orderings to output a summary file containing the values of delta for each progressive combination for each ordering. The output file will be named `small_fish_progu10_5_dashingsummary.csv`.
```
../lib/dandd progressive --dtree small_fish_5_dashing_dtree.pickle --norder 10 
```

This creates the internal file `sketchdb/small_fish_5_orderings.pickle` containing the orderings. It might be the case that these orderings are desired to be used for a parallel dataset (say the same samples stratified by variant type, as in the paper). Alternatively, this allows for reproduction of these initial orderings when additional orderings are added in a later run. In the command below another 20 orderings are stacked on top of the original 10 to produe a file called `small_fish_progu30_5_dashingsummary.csv`:

```
../lib/dandd progressive --dtree small_fish_5_dashing_dtree.pickle --norder 30 
```

### K-Independent Jaccard
If we wish to calculate the K-Independent Jaccard (KIJ) for our set, we use the `kij` command. If we also want the non-k independent pairwise jaccard values we use the `--jaccard` flag and provide a window of `--mink` and `--maxk` over which we would like our non-k independent Jaccards to be calculated.

```
../lib/dandd kij --dtree small_fish_5_dashing_dtree.pickle
```

### K-Sweep
In the paper we analyze what potential values of the delta ratio were considered in order to compare its growth to other measures as well as illustrate the selection of the argmax-k. In order to accomplish this, we used the `--ksweep` command to produce all potential values of delta for ks within a specific range. Again we provide a window of `--mink` and `--maxk`. Depending on which subcommand we provide we can sweep bare sketches with unions or every union performed in a progressive union. In the command below we do the latter:

```
../lib/dandd progressive --ksweep --dtree /home/jbonnie1/dandd/dev-dandD/example/small_fish_progu30_5_dashing_dtree.pickle
```
### Exactitude with KMC3
If you wish to count kmers rather than estimate with dashing, you would adjust your `tree` subcommand at the beginning to use the `--exact` argument. All files output after this change will replace the word "dashing" with "kmc" in their output filenames.