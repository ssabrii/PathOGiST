# PathOGiST

## Tasks
- [ ] Testing
  - [ ] Create complete set of unit test scripts for majority of functions
  - [ ] Add a continuous integration suite (such as [Travis CI](https://travis-ci.org/))
- [ ] Distribution
  - [x] Create [Conda](https://conda.io) package
  - [ ] Create [pip](https://pip.pypa.io/en/stable/) package
  - [ ] Create Galaxy version
- [ ] Lawyering up
  - [x] Add a license (e.g. MIT)
  - [ ] In each package file, add a header with copyright information
- [ ] Documentation

## Dependencies
- [`python3`](https://python.org)
- [`numpy`](https://numpy.org)
- [`scipy`](https://scipy.org)
- [`pandas`](https://pandas.pydata.org)
- [`matplotlib`](https://matplotlib.org)
- [`sklearn`](http://scikit-learn.org/stable/)
- [`pyyaml`](https://pyyaml.org/)
- [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) and the CPLEX Python API.
  Note: Do not install the Python CPLEX API by doing something like `pip install cplex`.
  This will install the trial version, which cannot handle large datasets.
  Instead, download CPLEX from the website, and run `python setup.py install` within the Python API folder.
  
## Installation
We recommend you create a conda environment for PathOGiST.
You can easily create one by running the following command in the top directory:
```bash
conda env create -f conda/environment.yaml
```
This will create a conda environment called `pathogist`, which you can activate by running
```bash
source activate pathogist
```
and deactivate with
```bash
source deactivate
```
When inside the `pathogist` conda environment, you can then simply run `/path/to/pathogist.py`.
Note that you will need to install CPLEX separately, as CPLEX is proprietary software.

## Subcommands

### Entire Pipeline (`all`)
This subcommand runs the PathOGiST pipeline from start to finish 
(i.e. distance matrix creation -> correlation clustering -> consensus clustering).

The main input file is a YAML configuration file, which you can create with the command
```bash
pathogist all [path to where you want your config] --new_config
```
The configuration file will look like so:
```bash
---
# PathOGiST configuration file.
# This configuration file is in YAML file format, release 1.2 (Third Edition).
# Google yaml for the specification.
# Do not remove the sections 'genotyping', 'distances', 'thresholds', 'all_constraints', or 'output'.
# Remove all key-value pairs in the sections 'genotyping' or 'distances' if you want them blank.
# Remove all list items in the section 'fine_clusterings' if you want that blank, too.
# Keys from 'genotyping' and 'distances' sections should not overlap.
genotyping:
    # Raw genotyping data from which to create distance matrices.
    # Accepted values are paths to text files containing paths to call files.
    # Currently only compatible with SNPs from snippy, MLSTs from MentaLiST, and CNVs from PRINCE.
    # Example:
    # MLST: /path/to/mlst_calls.txt
    SNP:
    MLST:
    CNV:
distances:
    # Paths to pre-constructed distance matrices in tsv format.
    # You can also specify SNP, MLST and CNV distance matrices here if you pre-constructed them,
    # but then they shouldn't appear in the section 'genotyping'.
    kmer:
    spoligotyping:
fine_clusterings:
    # The genotyping datatypes that are considered to be the "finest".
    - MLST
    - SNP
    - kWIP
thresholds:
    # Threshold values for performing correlation clustering on genotyping data types given above.
    # Every key appearing in the sections 'genotyping' and 'distances' should appear here with a value. 
    SNP:
    MLST:
    CNV:
    kmer:
    spoligotyping:
# Put the path for the final output consensus clustering (which will be in tsv format) here.
output:
# Use all constraints when performing correlation and consensus clustering
all_constraints: False
...
```
The inputs to the `genotyping` entries should be a file which contains absolute paths to your call files.
For example, `mlst_calls.txt` should look something like:
```bash
/path/to/SRR00001.calls
/path/to/SRR00002.calls
/path/to/SRR00003.calls
```
The output of PathOGiST is a TSV file containing the file consensus cluster assignment for each sample.

### Correlation Clustering (`correlation`)
The inputs to correlation clustering are:
* A distance matrix in the form of a TSV file
* A theshold cutoff value for the construction of the similarity matrix
The output is a TSV file containing the cluster assignments of the samples described by the distance matrix.

You can run correlation clustering with the following command:
```bash
pathogist.py correlation [distance matrix] [threshold] [output path]
```

### Consensus Clustering (`consensus`)
The input for consensus clustering is three files:
* A text file containing paths to distance matrices in `.tsv` format.
* A text file containing paths to clustering assignments in `.tsv` format.
* A text file containing the names of the clusterings which are 'finest'.
The output is a TSV file containing the cluster assignments of the samples which are common to all the input distance matrices.

You can run consensus clustering with the following command:
```bash
pathogist.py consensus [distances] [clusterings] [fine_clusterings] [output path]
```

Each line of the input files should correspond to a specific data type, e.g. SNPs, MLSTs, or CNVs.
Paths to distance matrices and cluster assignments should be prepended with the name of the clustering and an equal sign, i.e. `[name]=[path to file]`.
An example:

_Distances file_
```bash
SNP=/path/to/snp_dist
MLST=/path/to/mlst_dist
CNV=/path/to/cnv_dist
```
_Clusterings file_
```bash
SNP=/path/to/snp_clust
MLST=/path/to/mlst_clust
CNV=/path/to/cnv_clust
```
_Fine clusterings file_
```bash
SNP
```
