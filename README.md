# DandD  
A tool to estimate deltas for sequence sets and answer questions about relative contribution. DandD can be used to addresses how much new sequence is gained when a sequence collection grows and describe how much structural variation is discovered in each new human genome assembly allowing prediction of when discoveries will level off in the future. DandD uses a measure called $\delta$ (“delta”), developed initially for data compression and chiefly dependent on k-mer counts. DandD rapidly estimates $\delta$ using genomic sketches. Further information about the use of delta to address sequence similarity and how it is implemented in DandD can be found in the associated publication in iScience Volume 27, Issue 3, 15 March 2024, 109054.

 Jessica K. Bonnie, Omar Y. Ahmed, Ben Langmead. [**DandD: efficient measurement of sequence growth and similarity**]
 (https://www.sciencedirect.com/science/article/pii/S258900422400275X), [iScience](https://www.sciencedirect.com/journal/iscience) [10.1016/j.isci.2024.109054](https://doi.org/10.1016/j.isci.2024.109054)
 

## Installation  
```
git clone https://github.com/jessicabonnie/dandd  
cd dandd
chmod +x lib/dandd # make dandd executable
# NOTE: The following will alias dandd for the session. For this to persist this should be added to .bashrc
alias dandd=$(pwd)/lib/dandd 
```
### Dependencies
#### Dashing
**Note**: The version of *Dashing* currently available on bioconda is **not** recent enough. Please obtain the latest binary release for your system from [https://github.com/dnbaker/dashing-binaries](https://github.com/dnbaker/dashing-binaries) and move it to somewhere in your `$PATH`. Alternatively, you can install the latest version from source by cloning it from github
: [https://github.com/dnbaker/dashing](https://github.com/dnbaker/dashing):  
```
git clone --recursive https://github.com/dnbaker/dashing  
cd dashing && make dashing
```
#### KMC3
If you wish to use the `--exact` command to retrieve the *actual* k-mer counts rather than using Dashing to estimate them, you will need to install `kmc` (with `kmc_tools`) and make sure that it is in your `$PATH`. You can install it through bioconda with the following command:
```
conda install -c bioconda kmc
```
#### Parallel
This version of DandD uses `parallel` to simultaneously create sketches across different values of `k`. You will need to make sure that your system has this by loading the module on your server or installing it on your computer.  
```
conda install -c bioconda parallel
```

## File Output 
The tool will write summary files to the current working directory or to the output directory provided. If a sketch directory is not provided it will be sought or created within the output directory. The sketch directory will contain (1) files mapping the sketch/database names to the component inputs (2) subdirectories of sketches/dbs arranged by number of genomes and k-mer length.

**outputdir /**  
* |_ deltatree : pickle to be used internally during later calls to DandD 
* |_ deltas and summary tables : user readable tables named using tag/labels provided   
* |_ sketchdb.txt tracking file : user readable table tracking locations and inputs for sketches and database files relevant to particular experiment (prefixed accordingly)  
* **sketchdb /**  
   * |_ dandd_fastahex, dandd_sketchinfo, dandd_cardinalities : internal pickles tracking name assignments and cardinalities.  
   * **ngen[#] /** : one directory for each number of genomes/assemblies/fastas unioned (starting at 1 for single fasta sketches/databases)     
     - **k[#] /** : one directory for each k explored  
       * |_ dashing-sketch  
       * |_ kmc-database files 


## Commands

### Programmatic
There are a couple of commands options that pertain to the program in general. These flags must be provided before the `tree`, `kij`, etc. subcommands.  
`--debug`: this flag will cause DandD to print the commands for each call to Dashing or KMC so that you can run them directly to see what kind of error they are returning.  
`--fast`: Avoid writing out intermediate copies of the reference files. (Meaning that cardinalities and pathing information will not be saved if DandD is terminated prematurely.)  
`--safe`: Double check all sketch hashes to make sure they match the sums of the fasta hashes. This is time consuming, but is recommended if you have any suspicion that the fasta files have changed.  
`--verbose`: Un-supress messages from DandD about progress. Unfortunately, there is no way currently to supress output from Dashing and KMC regardless of this flag.  
### tree
The `tree` command is the first step in performing any analysis. It creates a pickle with the structure DandD will expect to receive for the other commands in order to answer questions downstream. It **requires** either a file containing full filepaths to fastas (`--fastas <path>`, `-f <path>`) or a data directory (`--datadir <path>`) containing all fastas of interest.  
`--outdir <path>` or `-o <path>` : output directory (default current working directory)  
`--sketchdir <sketchdb path>` or `-c <sketchdb path>` : path to sketch directory. If the sketches for `tree` already exist in a `sketchdb` directory other than inside `outdir`, the `sketchdir` can be provided in order to take maximum advantage of previous work in terms of time/space. This is especially important for the sketches in `sketchdb/ngen1` which take the most time and space when using Dashing.  
Additionally:  
`-s <string>` or `--tag <string>` : string used in file name prefixes to name the non-sketch output files for a particular experiment. If no tag is provided, `"dandd"` will be used by default.  
`-k <int>` : set the starting value in search for optimal `k` for delta. As demonstrated in the DandD paper, the optimal k-mer length for delta changes based on characteristics of the data. Therefore, it is beneficial to adjust the starting value of k based on prior information to avoid searching excessively high or low values.  
`--exact` : use KMC for k-mer counting (instead of estimation through Dashing sketching) to calculate delta.   
`--non-canon` : disable canonicalization. By default, DandD will canonicalize the k-mers (i.e. treat reverse complements as equivalent to the original k-mer string).  
`--nchildren <int>` or `-n <int>` : maximum number of children for each union. By default, DandD will build trees with a single layer of child nodes and one single union of all inputs (nicknamed 'spiders'). With `--nchildren` the tree can be built by calculating delta for a series of smaller unions by indicating the maximum number of children for each union.  
`--registers <int>` : adjust number of registers used during sketching. For information about the implementation and trade-offs of adjusting this value, see Dashing's documentation.  
`--lowmem` : delete union sketches/databases after storing their cardinalities and use stored values even when the sketch is missing. The path information etc. for the sketch will remain stored in the internal files such as dandd_fastahex and dandd_sketchinfo as well as the sketchdb.txt for the experiment.  **Note:** It is not recommended to use this command the first time the databases are initialized.  
#### Result Files  
Result files are named using a **prefix** composed of: the provided `tag` string, the total number of fasta files used to produce them, and the `tool` used (kmc or dashing).  
**<prefix>_deltas.csv**  
The `tree` command outputs the deltas for the component fastas and their full union as well as any deltas of intermediate unions if the `--step <int>` argument is used. (In addition to producing the pickle file to be used by other subcommands.) The field order is subject to change but columns contain the following:  
*  ngen - the number of fasta files in the union
*  title - the basenames of each fasta separated by underscores
*  sketchloc - the path to the sketch whose cardinality was used to find delta
*  fastas - the full paths to the component fastas separated by "|"
*  card - the cardinality (count or estimated count) of unique k-mers of length k
*  delta - the value of delta for the union
### progressive

The `progressive` command performs a series of cumulative (or progressive) unions over a number (`--norder`) of orderings of the input fastas. It **requires** `--dtree / -d <tree.pickle>` where `<tree.pickle` is an output of the `tree` command. It will automatically use the same values for `--lowmem`, `--step`, `--sketchdir`, `--non-canon`, `--registers` embedded in that delta-tree output. If only a subset of the original fastas in the tree are of interest, a subset list can be provided using `--fastas / -f <filepath>`. If orderings have previously been produced for this tree, those orderings will be repurposed, with additional orderings generated if `--norder` is greater than the number of orderings already generated.  If an ordering file is to be shared across experiments/tags, a preexisting ordering file can be used via `--orderings / -r <ordering pickle>`.  If only one deterministic ordering is desired provide it via the fasta file and use `--norder 1`.
`--ksweep`: populate all of the possible delta values across a range of `k` values for each sketch/database combination in the progressive union. `--mink <int>` and `--maxk <int>` are used to bound that range with default [2,32]. **Note** For reasons internal to Dashing it must be the case that maxk<= 32 for estimation. If you want to sweep higher `k`s, you must use KMC via the `--exact` flag. If `--ksweep` is set, the progressive output will not seek or produce the delta values for each union, but instead provide the *possible delta* values (output field: delta_pos) for each union at each specified k in the range.

#### Result Files
Result files are named using a **prefix** composed of: the provided `tag` string, the total number of fasta files used to produce them, progu<the number of orderings (from `--norder` or zero if entirety of default ordering file is used)>, and the `tool` used (kmc or dashing).  
**<prefix>.csv**  
If `--ksweep` is **not** included, this file will contain the `progressive` outputs of the deltas for the component fastas for each union in each ordering permutation. If `--ksweep` is provided but the argmax-k is not in [mink, maxk], it will be nonsense. The field order is subject to change but columns contain the following:  
*  ordering - the index of the ordering (one-based numbering)  
*  ngen - the number of fasta files in the union; also the index *within* the ordering  
*  title - the basenames of each fasta separated by underscores
*  sketchloc - the path to the sketch whose cardinality was used to find delta
*  fastas - the full paths to the component fastas separated by "|"
*  card - the cardinality (count or estimated count) of unique k-mers of length k
*  delta - the value of delta for the union

**<prefix>_<tool>summary.csv**  
This file contains all of the possible delta values across all *k*s scanned (whether during the search for delta or via `--ksweep`) during every step of every ordering. The field order is subject to change but columns contain the following:  

*  ordering - the index of the ordering (one-based numbering)  
*  ngen - the number of fasta files in the union; also the index *within* the ordering  
*  kval - the k-value (k-mer length) that was used to find the cardinality and the delta-pos  
*  title - the basenames of each fasta in the union separated by underscores
*  command - the bash command used to produce the sketch in case you wish to run the individual command (n.b the component ngen1 sketches will need  to still be at the indicated location)  
*  card - the cardinality (count or estimated count) of unique k-mers of length kval  
*  delta_pos - the possible value of delta for this union at this kval (i.e. card/k)  
### kij
The `kij` command produces pairwise deltas across all of the input fastas in the `tree.pickle` in order to compute K-Independent Jaccards. It also calculates values for pairwise Jaccard within the range provided using `--mink` and `--maxk`. To write out those pairwise values use `--jaccard`. It **requires** `--dtree / -d <tree.pickle>` where `<tree.pickle` is an output of the `tree` command.  

#### Result Files
Result files are named using a **prefix** composed of: the provided `tag` string, the total number of fasta files used to produce them, progu<the number of orderings (from `--norder` or zero if entirety of default ordering file is used)>, and the `tool` used (kmc or dashing).  
**<prefix>.kij.csv**  
ABk,Btitle,ABdelta,Atitle,Bdelta,A,Bk,Adelta,Ak,KIJ,B
This file contains the values for the k-Independent Jaccard similarity metric.  
*  <A,B> - file path of FASTA <A,B>  
*  <A,B>title - stripped sample id from FASTA <A,B>  
*  <A,B>delta - value of delta for <A,B> alone  
*  <A,B>k - argmax-k used to calculate <A,B>delta
*  ABdelta - value of delta for union of A and B
*  ABk - argmax-k used to calculate ABdelta
*  KIJ - value of k-Independent Jaccard similarity metric

**<prefix>.j.csv**  
This file contains the values for the standard Jaccard similarity metric if the `--jaccard` flag and `mink`/`maxk` values are provided.  

**<prefix>_AFtuples.pickle**  
This file contains a pickle that can be passed to unsupported functions in `helpers/afproject.py`.
