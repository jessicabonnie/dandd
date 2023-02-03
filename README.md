# DandD  
Tool to estimate deltas for sequence sets and answer questions about relative contribution

## Installation
**Note**: The version of *Dashing* currently available on bioconda is **not** recent enough. Please install the latest version from github
: [https://github.com/dnbaker/dashing](https://github.com/dnbaker/dashing)


## File Output 
The tool will write summary files to the current working directory or to the output directory provided. If a sketch directory is not provided it will be sought or created within the output directory. The sketch directory will contain (1) files mapping the sketch/db names to the component inputs (2) subdirectories of sketches/dbs arranged by number of genomes and k-mer length.

outputdir/  
*  |_ tree pickle  
*  |_ summary files (named using tag/labels provided)  
*  |_ sketchdb/  
   * |_ dandd_fastahex.pickle  
   * |_ dandd_sketchinfo.pickle  
   * |_ dandd_cardinalities.pickle  
   * |_ ngen1/  
     * |_ k{NUM}  
   * |_ ngen{NUM}  
