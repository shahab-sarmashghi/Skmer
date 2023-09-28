[![Anaconda-Server Badge](https://anaconda.org/bioconda/skmer/badges/version.svg)](https://anaconda.org/bioconda/skmer)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skmer/badges/downloads.svg)](https://anaconda.org/bioconda/skmer)

# Skmer
Skmer is a fast tool for estimating distances between genomes from low-coverage sequencing reads (genome-skims), without needing any assembly or alignment step. The paper where we have described the methods and tested Skmer on simulated short reads and SRA's from previous sequencing experiments is available online (open access):
  - [Sarmashghi, S., Bohmann, K., P. Gilbert, M. T., Bafna, V., & Mirarab, S. (2019). Skmer: assembly-free and alignment-free sample identification using genome skims. Genome Biology, 20(1), 34. https://doi.org/10.1186/s13059-019-1632-4][1]

And the paper where we have described the procedure for estimating branch support for phylogenies generated using Skmer is here:

 - [Eleonora Rachtman, Shahab Sarmashghi, Vineet Bafna, Siavash Mirarab, Quantifying the uncertainty of assembly-free genome-wide distance estimates and phylogenetic relationships using subsampling, Cell Systems, Volume 13, Issue 10, 2022, Pages 817-829.e3, ISSN 2405-4712, https://doi.org/10.1016/j.cels.2022.06.007.][8]

Skmer is a command-line tool implemented in python. It runs [Jellyfish][2] and [Mash][3] internally to efficiently compute k-mer profile of genome-skims and their intersection, and estimates the genomic distances by correcting for the effect of low coverage and sequencing error. Skmer also depends on [seqtk][5] for some FASTQ/A processings. After following the installation steps and going over the usage guide, check out these awesome [tutorials][9] by Siavash Mirarab to learn more about using Skmer and other tools we have developed for genome skimming.

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
Skmer has five sub-commands:

### reference
Gets the path to a directory of FASTQ/FASTA files (one uncompressed *.fastq/.fq/.fa/.fna/.fasta* file per each sample) and creates a reference library containing the estimates of sequencing parameters as well as the Mash sketch for each sample. If you have paried-end sequencing reads, we suggest using tools such as [BBMerge][10] to merge overlapping read pairs. If the input is an assembled sequence (determined by the length of sequences) the correction for low coverage and sequencing error is not applied to that sample. All corrected pairwise genomic distances are then estimated and written to a file. For a test run, change the directory to `data` under your Skmer installation directory, and run
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
Gets the path to a directory of FASTQ/FASTA files (one uncompressed *.fastq/.fq/.fa/.fna/.fasta* file per each sample) and performs subsampling procedure. Function will create `subsample` directory containing replicate subfolders `rep0`, `rep1` etc.
```
skmer subsample ref_dir
```
A number of additional paramters can be specified. `-S` allows to provide custom seed that will be used to generate a list of seeds for each subreplicate (default is 42). With option `-sub` the user can define directory of output for subsample replicates (default is `working_directory/subsample`). `-b` option can be used to indicate subreplicate count (by default value is set to 100). `-i` allows to specify index of the first replicate (default is 0). 
```
skmer subsample -b 100 ref_dir -s 100000 -S 42 -p 24 -t -i 0
```
To see the complete list of inputs and options, run `skmer subsample -h`.

#### Helpful tips:
- **Job parallelization:** Combinations of `-b` and `-i` allows a flexible job parallelization. For example, to test large dataset user can run subsampling in chunks by specifying `-b 10 -i 0 -S 14500` (generates 10 subreplicates starting at index 0 such as first repository is rep0 and others are rep1, rep2 ... rep9), `-b 10 -i 10 -S 13800` (generates 10 replicates starting at index 10 such as subrepicates are rep10, rep11 ... rep19). Here we note that since internally Skmer uses default seed 42 when subsampling job is split **_variable seed is needed to be specified_** otherwise replicates will come out the same.
- **Sample uniformity:** At the moment subsample only works for cases where all samples are either assemblies or sequencing reads. If there is a need to run a combination the user can simulate sequencing reads from assemblies. We recommnd tool such as [ART][6] for this purpose.

### correct
Performs correction of subsampled distance matrices obtained for reference genome-skims or assemblies. Since distance matrices are precomputed this step is fast. 
Output is this command is a set of corrected distance matrices for main estimate and subreplicates. 

Input main distance matrix can have any filename. After correction new matrix file will be created with the name appended with the suffix `_cor_`. We note that main distance matrix remains unchanged and correction only involves rounding of the values up to 12 significant digits to ensure that output is compatible with downstream tools like FastMe.

Correction algorithm looks for `dimtrx_rep.txt` file in each subreplicate directory. For every subreplicate matrices with both types of correction are generated. Corrected distance matrices are appended with suffixes `_cor` and `_cor_cons` which corresponds to main and consensus correction respectively.
```
skmer correct -main jc-dist-mat -sub subsample_dir
```
`-main` option takes as an input distance matrix file for main estimate before subsampling. This should be computed using standard `reference` command. `-sub` is used to specify location of `subsample` directory. These options have no default settings.

<br/>

Workflow for computing _k-mer_-based trees with branch support
-----

### 1. Getting Skmer distance matrices
We suggest the following workflow to obtain Skmer distance matrices for sequencing reads or assemblies.

* **To obtain main estimate distance matrix before subsampling:**
```
skmer reference ref_dir -s 100000 -S 42 -p 24 -t -o dimtrx_main
```
* **To generate subreplicates:**  
```
skmer subsample -b 100 ref_dir -s 100000 -S 42 -p 24 -t -i 0
```
* **To correct estimates:**
```
skmer correct -main path_to_file/dimtrx_main.txt -sub path_to_directory/subsample
```

### 2. Reformatting distance matrices
In order to be compatible with downstream software, standard square distance matrices should be converted into [PHYLIP][7] format.

Example of distance matrix in standard square form
```
sample     Alpha Beta  Gamma Delta Epsilon
Alpha      0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 3.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
```

Example of distance matrix in PHYLIP format
```
    5
Alpha      0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 3.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
```
One way to perform the formatting is to run our custom [script](https://github.com/noraracht/Skmer/blob/master/data/tsv_to_phymat.sh) that is available in data folder:
```
bash tsv_to_phymat.sh dimtrx_original.txt dimtrx_reformatted.txt
```

### 3. Phylogeny inference
[FastME](http://www.atgc-montpellier.fr/fastme/) can be used to infer phylogenies from reformatted distance matrices. With default parameters the program can be launched: 
```
fastme -i input_data_file -o output_tree_file
```
In our case input_data_file contains reformatted distance matrix and output_tree_file contains computed backbone tree.

### 4. Concatenation of trees
To concatenate output phylogenies for subreplicates the user can run the following commands:

**For main correction**
```
cat subsample/rep*/dimtrx_rep_cor.txt.tre > bootstrap.All
```
**For consensus**
```
cat subsample/rep*/dimtrx_rep_cor_cons.txt.tre > bootstrap.All_consensus
```

### 5. Generation of final tree with estimated support
To generate phylogeny with support values we suggest to use [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/).

To compute phylogeny using **main correction** the user can draw bipartition information on a tree provided with -t (e.g., main Skmer estimate) based on multiple trees (e.g., from subsampled replicates) in a file specified by -z using command:
```
raxmlHPC -f b -m GTRCAT -z bootstrap.All -t dimtrx_main_cor_.txt.tre -n BS_TREE_MAIN
```
To compute extended majority rule **consensus** tree with "-J MRE" out of the subsampled replicates the user can run:
```
raxmlHPC -J MRE -z bootstrap.All_consensus -p 4424 -m GTRCAT -n BS_TREE_CONS
```

Note: some visualization tools, for instance [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), are not compatible with the consensus output. If there is a need to make it work, phylogeny file can be reformatted using:
```
sed -E 's/([:][0-9]+[.][0-9]+)[[]([0-9]+)[]]/\2\1/g' RAxML_MajorityRuleExtendedConsensusTree.BS_TREE_CONS > RAxML_MajorityRuleExtendedConsensusTree.BS_TREE_CONS_fixed
```

Tutorials
---------
Check out the following [tutorials][9] for more on how to use Skmer and other tools we have developed for different genome skimming applications.

[1]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1632-4
[2]: http://www.genome.umd.edu/jellyfish.html
[3]: http://mash.readthedocs.io/en/latest/
[4]: https://conda.io/miniconda.html
[5]: https://github.com/lh3/seqtk
[6]: https://manpages.debian.org/testing/art-nextgen-simulation-tools/art_illumina.1.en.html
[7]: https://evolution.genetics.washington.edu/phylip/doc/distance.html
[8]: https://www.sciencedirect.com/science/article/pii/S2405471222002770
[9]: https://github.com/smirarab/tutorials
[10]: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/

