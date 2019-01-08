# Skmer
Skmer is a fast tool for estimating distances between genomes from low-coverage sequencing reads (genome-skims), without needing any assembly or alignment step. The preprint where we have described the methods and tested Skmer on simulated short reads is available on bioRxiv:
  - [Sarmashghi, S., Bohmann, K., Gilbert, M. T. P., Bafna, V., & Mirarab, S. (2018). Assembly-free and alignment-free sample identification using genome skims. bioRxiv. https://doi.org/10.1101/230409][1]

Skmer is a command-line tool implemented in python. It runs [Jellyfish][2] and [Mash][3] internally to efficiently compute k-mer profile of genome-skims and their intersection, and estimates the genomic distances by correcting for the effect of low coverage and sequencing error. 

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
2. Install [Jellyfish][2] (v2.2.6 or later) and [Mash][3] (v1.1 or later), and add the path to their binary to the system path (so you can, e.g., run `jellyfish --version` and `mash --version` successfully in the terminal). 
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
Skmer has two sub-commands:

### reference
Gets the path to a directory of FASTQ files (one uncompressed *.fastq* genome-skim per each species), and creates a reference library containing the estimates of sequencing parameters as well as the Mash sketch for each genome-skim. All corrected pairwise genomic distances are estimated and written to a file. For a test run, change the directory to `data` under your Skmer installation directory, and run
```
skmer reference ref_dir 
```
The genome-skims in `ref_dir` directory are processed, and a reference `library` is created in the working directory. You can specify a custom name (and so its path) for your library using `-l` option
```
skmer reference ref_dir -l custom_library_name
```
Default k-mer size is set to `31` which is the maximum length allowed by Mash, and can be changed using `-k` option. We do not recommend using k-mers smaller than ~`21`, as k-mers without any shared evolutionary history start to seem similar just out of random chance. The sketch size can also be changed using `-s` option from its default value `10000000`. Decreasing the sketch size will reduce the size of library on disk, but also compromises the accuracy of distance estimation. The corrected pairwise distances are estimated and written to the file `ref-dist-mat.txt` in the working directory by default. The output prefix (and so its path) can be changed using `-o` option
```
skmer reference ref_dir -o output_prefix
```
If distances are going to be used to build phylogenies, it is recommended to use `-t` flag. In this case, the estimated distances are transformed to the phylogenetic distances using the Jukes-Cantor model of substitution. Run `skmer reference -h` for the complete list of arguments and options.  

### query
Processes a query genome-skim, and outputs the sorted list of reference species based on their distance to the query. Optionally, the query can be added to the reference library. To test its function, assuming that you have already run `reference` and compiled the reference `library`, in `data` directory run
```
skmer query qry.fastq library
```
The sorted list of reference species and their distances from the query is written to `dist-qry.txt`. You can change the output prefix from `dist` to something else using `-o` option
```
skmer query qry.fastq library -o output_prefix
```
If you want to add the processed query to the reference library and include it as a reference for future comparisons, use `-a` flag. To see the complete list of inputs and options, run `skmer query -h`.

[1]: https://www.biorxiv.org/content/early/2018/04/02/230409
[2]: http://www.genome.umd.edu/jellyfish.html
[3]: http://mash.readthedocs.io/en/latest/
[4]: https://conda.io/miniconda.html
