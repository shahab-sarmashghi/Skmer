# Skmer (Under construction)
Skmer is a fast tool for estimating distances between genomes from low-coverage sequencing reads (genome-skims), whithout needing any assembly or alignment step. The preprint where we have described the methods and tested Skmer on simulated short reads is available on bioRxiv:
  - [Sarmashghi, S., Bohmann, K., Gilbert, M.T.P., Bafna, V., Mirarab, S.: Assembly-free and alignment-free sample identification using genome skims. bioRxiv (dec 2017) 230409][1]

Skmer is a command-line tool implemented in python. It runs [Jellyfish][2] and [Mash][3] internally to efficiently compute k-mer profile of genome-skims and their intersection, and estimates the genomic distances by correcting for the effect of low coverge and sequencing error. 

Installation
------------
**On 64-bit Linux and Mac OSX**, you can install Skmer from bioconda channel using conda package manager. 
1. Install [Miniconda][4] (you can skip this if you already have either of Miniconda or Anaconda installed). 
2. Add the bioconda channel by running the following commands in your terminal (order matters):
```
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
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

Usage
------------
Skmer has two sub-commands:
### reference
Gets the path to a directory of FASTQ files (one uncompressed *.fastq* genome-skim per each species), and creates a reference library containing the estimates of sequencing parameters as well as the Mash sketch for each genome-skim. All corrected pairwise genomic distances are estimated and written to a file. Run `skmer reference -h` for detailed information about the arguments and options.
### query
Processes a query genome-skim, and outputs the sorted list of reference species based on their distance to the query. Optionally, the query can be added to the reference library. Run `skmer query -h` for a complete list of inputs and options.

[1]: https://www.biorxiv.org/content/early/2017/12/08/230409
[2]: http://www.genome.umd.edu/jellyfish.html
[3]: http://mash.readthedocs.io/en/latest/
[4]: https://conda.io/miniconda.html
