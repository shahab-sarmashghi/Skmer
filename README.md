[![Anaconda-Server Badge](https://anaconda.org/bioconda/skmer/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skmer/badges/downloads.svg)](https://anaconda.org/bioconda/skmer)

# Skmer
Skmer is a fast tool for estimating distances between genomes from low-coverage sequencing reads (genome-skims), without needing any assembly or alignment step. The paper where we have described the methods and tested Skmer on simulated short reads and SRA's from previous sequencing experiments is available online (open access):
  - [Sarmashghi, S., Bohmann, K., P. Gilbert, M. T., Bafna, V., & Mirarab, S. (2019). Skmer: assembly-free and alignment-free sample identification using genome skims. Genome Biology, 20(1), 34. https://doi.org/10.1186/s13059-019-1632-4][1]

The paper where we have described **_procedure for estimating branch support for phylogenies generated using Skmer_** will be available online shortly. We are working on integrating changes into the main Skmer branch.

Skmer is a command-line tool implemented in python. It runs [Jellyfish][2] and [Mash][3] internally to efficiently compute k-mer profile of genome-skims and their intersection, and estimates the genomic distances by correcting for the effect of low coverage and sequencing error. Skmer also depends on [seqtk][5] for some FASTQ/A processings. 

Installation
------------
**On 64-bit Linux and Mac OSX**, you can install Skmer from bioconda channel using conda package manager. 
1. Install [Miniconda][4] (you can skip this if you already have either of Miniconda or Anaconda installed). 
2. Add the bioconda channel by running the following commands in your terminal (order matters):
```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
```
3. Run the following command to install Skmer (and all dependencies) 
```
    conda install skmer
```

**Alternatively, and for all other OS and architectures**, you can download the github repository and install Skmer using the setup script. 
1. You need to have python 2.7 or later installed
2. Install [Jellyfish][2] (v2.2.6 or later), [Mash][3] (v2.3 or later), and [seqtk][5] (v1.3), and add the path to
 their binary to the system path (so you can, e.g., run `jellyfish --version`, `mash --version`, and `seqtk` successfully in the terminal). 
3. Clone the github repository by running (or you can download the repo)
```
    git clone https://github.com/shahab-sarmashghi/Skmer.git
```
4. Change to the Skmer directory and run
```
    python setup.py install
```

Using Skmer
------------
Skmer has three sub-commands:

### reference
Gets the path to a directory of FASTQ/FASTA files (one uncompressed *.fastq/.fq/.fa/.fna/.fasta* file per each sample) and creates a reference library containing the estimates of sequencing parameters as well as the Mash sketch for each sample. If the input is an assembled sequence (determined by the length of sequences) the correction for low coverage and sequencing error is not applied to that sample. All corrected pairwise genomic distances are then estimated and written to a file. For a test run, change the directory to `data` under your Skmer installation directory, and run
```
skmer reference ref_dir -p 4
```
The genome-skims and assemblies in `ref_dir` directory are processed (using 4 cores in parallel), and a reference `library` is created in the working directory. You can specify a custom name (and so its path) for your library using `-l` option
```
skmer reference ref_dir -l custom_library_name
```
Default k-mer size is set to `31` which is the maximum length allowed by Mash, and can be changed using `-k` option. We do not recommend using k-mers smaller than ~`21`, as k-mers without any shared evolutionary history start to seem similar just out of random chance. The sketch size can also be changed using `-s` option from its default value `10000000`. Decreasing the sketch size will reduce the size of library on disk, but also compromises the accuracy of distance estimation. The corrected pairwise distances are estimated and written to the file `ref-dist-mat.txt` in the working directory by default. The output prefix can be changed using `-o` option
```
skmer reference ref_dir -o output_prefix
```
If distances are going to be used to build phylogenies, it is recommended to use `-t` flag. In this case, the estimated distances are transformed to the phylogenetic distances using the Jukes-Cantor model of substitution. Run `skmer reference -h` for the complete list of arguments and options.  


### distance
Computes all pairwise distances for a processed library. The main usage is to compute distances when combining already processed libraries, otherwise `reference` command outputs distances as well when the input files are processed. For example, in `data` directory, assuming that you have already run `reference` and compiled the reference `library`, try
```
skmer distance library -t -o jc-dist-mat
```
The distances between all samples in the `library` are computed, and after applying Jukes-Cantor transformation, the mutation rates are written to `jc-dist-mat.txt`, which can be later used to build a phylogeny based on distances. To see the help message, run `skmer distance -h`.

### query
Processes a query genome-skim or assembly, and outputs the sorted list of reference samples based on their distance to the query. Optionally, the query can be added to the reference library. To test its function, assuming that you have already run `reference` and compiled the reference `library`, in `data` directory run
```
skmer query qry.fastq library
```
The sorted list of reference species and their distances from the query is written to `dist-qry.txt`. You can change the output prefix from `dist` to something else using `-o` option
```
skmer query qry.fastq library -o output_prefix
```
If you want to add the processed query to the reference library and include it as a reference for future comparisons, use `-a` flag. To see the complete list of inputs and options, run `skmer query -h`.

### subsample
Gets the path to a directory of FASTQ/FASTA files (one uncompressed *.fastq/.fq/.fa/.fna/.fasta* file per each sample) and performs subsampling procedure. Function will create `subsample` folder containing replicate subfolders `rep0`, `rep1` etc.
```
skmer subsample ref_dir
```
A number of additional paramters can be specified. `-b` option can be used to indicate subreplicate count (by default value is set to 100). `-i` allows to specify index of the first replicate (default is 0). Combinations of `-b` and `-i` should allow for a more flexible job parallelization. `-S` allows to provide custom seed that will be used to generate a list of seeds for each subreplicate (default is 42). With option `-sub` the user can define directory of output for subsample replicates (default is `working_directory/subsample`)
```
skmer subsample -b 100 ref_dir -s 100000 -S 42 -p 24 -t -i 0
```
To see the complete list of inputs and options, run `skmer subsample -h`.

### correct
Performs correction of subsampled distance matrices obtained for reference
genome-skims or assemblies
```
skmer correct ref_dir/subsample di_matrx
```
To see the complete list of inputs and options, run `skmer correct -h`.

[1]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1632-4
[2]: http://www.genome.umd.edu/jellyfish.html
[3]: http://mash.readthedocs.io/en/latest/
[4]: https://conda.io/miniconda.html
[5]: https://github.com/lh3/seqtk
