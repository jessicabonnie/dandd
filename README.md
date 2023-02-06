# DandD  
Tool to estimate deltas for sequence sets and answer questions about relative contribution

## Installation
### Dependencies
#### Dashing
**Note**: The version of *Dashing* currently available on bioconda is **not** recent enough. Please install the latest version from github
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
The tool will write summary files to the current working directory or to the output directory provided. If a sketch directory is not provided it will be sought or created within the output directory. The sketch directory will contain (1) files mapping the sketch/db names to the component inputs (2) subdirectories of sketches/dbs arranged by number of genomes and k-mer length.

outputdir/  
* |_ tree pickle  
* |_ summary files (named using tag/labels provided)  
* |_ sketchdb/  
   * |_ dandd_fastahex.pickle  
   * |_ dandd_sketchinfo.pickle  
   * |_ dandd_cardinalities.pickle  
   * |_ ngen1/  
     - |_ k[num]  
  * |_ ngen[num] 


## Commands

### Programmatic
There are a couple of commands options that pertain to the program in general. These flags must be provided before the `tree`, `kij`, etc. subcommands.  
`--debug`: this flag will cause DandD to print the commands for each call to Dashing or KMC so that you can run them directly to see what kind of error they are returning.  
`--fast`: Avoid writing out intermediate copies of the reference files. (Meaning that cardinalities and pathing information will not be saved if DandD is terminated prematurely.)  
`--safe`: Double check all sketch hashes to make sure they match the sums of the fasta hashes. This is time consuming, but is recommended if you have any suspicion that the fasta files have changed.  
### tree
The `tree` command is the first step in performing any analysis. It creates a pickle with the structure DandD will expect to receive for the other commands in order to answer questions downstream. It requires either a file containing full filepaths to fastas (`--fastas`) or a data directory (`--datadir`) containing all fastas of interest.  
If the sketches for `tree` already exist in a `sketchdb` directory other than inside `--outdir` (default current working directory), the `--sketchdir` can be provided in order to take maximum advantage of previous work in terms of time/space.  

### progressive

The `progressive` command performs a series of cumulative (or progressive) unions over a number (`--norder`) of orderings of the input fastas. It takes a `tree.pickle` object output from the `tree` command. If only a subset of the original fastas in the tree are of interest, a subset list can be provided. If orderings have previously been produced for this tree, those orderings will be repurposed, with additional orderings generated if `--norder` is greater than the number of orderings already generated.  

### kij
The `kij` command produces pairwise deltas across all of the input fastas in the `tree.pickle` in order to compute K-Independent Jaccards. It also calculates values for pairwise Jaccard within the range provided using `--mink` and `--maxk`.  


### info
The `info` command is used to obtain `ksweep` and `summary` information.  
