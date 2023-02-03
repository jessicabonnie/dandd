# Example: Fish Mitochondria

**Note**: this assumes that `dashing` is in your `$PATH`

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
../lib/dandd_cmd.py tree --datadir ./assembled-fish_mito --tag fish-mito

```

## SubTree
Suppose you are only interested in a subset of the fastas for a particular experiment? You can do that too! DandD will use the sketches already produced to create this smaller tree.

```
realpath ./assembled-fish_mito/* | head > fish_subset.txt
../lib/dandd_cmd.py tree --fastas fish_subset.txt --tag small_fish -k 12

```

## Progressive Union

```

```
