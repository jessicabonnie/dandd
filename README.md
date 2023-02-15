# DandD  
Tool to estimate deltas for sequence sets and answer questions about relative contribution

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
`--sketchdir <sketchdb path>` or `-c <sketchdb path>` : path to sketchdirectory. If the sketches for `tree` already exist in a `sketchdb` directory other than inside `outdir`, the `sketchdir` can be provided in order to take maximum advantage of previous work in terms of time/space.  
Additionally:  
`-s <string>` or `--tag <string>` : string used in file name prefixes to name the non-sketch output files for a particular experiment. If no tag is provided, `"dandd"` will be used by default.  
`-k <int>` : set the starting value in search for optimal `k` for delta. As demonstrated in the DandD paper, the optimal k-mer length for delta changes based on characteristics of the data. Therefore, it is beneficial to adjust the starting value of k based on prior information to avoid searching excessively high or low values.  
`--exact` : use KMC for k-mer counting (instead of estimation through Dashing sketching) to calculate delta.   
`--non-canon` : disable canonicalization. By default, DandD will canonicalize the k-mers (i.e. treat reverse complements as equivalent to the original k-mer string).  
`--nchildren <int>` or `-n <int>` : maximum number of children for each union. By default, DandD will build trees with a single layer of child nodes and one single union of all inputs (nicknamed 'spiders'). With `--nchildren` the tree can be built by calculating delta for a series of smaller unions by indicating the maximum number of children for each union.  
`--registers <int>` : adjust number of registers used during sketching. For information about the implementation and trade-offs of adjusting this value, see Dashing's documentation.  


### progressive

The `progressive` command performs a series of cumulative (or progressive) unions over a number (`--norder`) of orderings of the input fastas. It **requires** `--dtree / -d <tree.pickle>` where `<tree.pickle` is an output of the `tree` command. If only a subset of the original fastas in the tree are of interest, a subset list can be provided using `--fastas / -f <filepath>`. If orderings have previously been produced for this tree, those orderings will be repurposed, with additional orderings generated if `--norder` is greater than the number of orderings already generated.  If an ordering file is to be shared across experiments/tags, a preexisting ordering file can be used via `--orderings / -r <ordering pickle>`.  

### kij
The `kij` command produces pairwise deltas across all of the input fastas in the `tree.pickle` in order to compute K-Independent Jaccards. It also calculates values for pairwise Jaccard within the range provided using `--mink` and `--maxk`. To write out those pairwise values use `--jaccard`. It **requires** `--dtree / -d <tree.pickle>` where `<tree.pickle` is an output of the `tree` command.  

### info
The `info` command is used to obtain `ksweep` and `summary` information. It **requires** `--dtree / -d <tree.pickle>` where `<tree.pickle` is an output of the `tree` command.  
