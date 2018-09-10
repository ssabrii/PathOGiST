# PathOGiST

## Tasks
- [ ] Testing
  - [ ] Create complete set of unit test scripts for majority of functions
  - [ ] Add a continuous integration suite (such as [Travis CI](https://travis-ci.org/))
- [ ] Distribution
  - [ ] Create [Conda](https://conda.io) package
  - [ ] Create [pip](https://pip.pypa.io/en/stable/) package
- [ ] Layering up
  - [ ] Add a license (e.g. MIT)
  - [ ] In each package file, add a header with copyright information
  
## Dependencies
- [`python3`](https://python.org)
- [`numpy`](https://numpy.org)
- [`scipy`](https://scipy.org)
- [`pandas`](https://pandas.pydata.org)
- [`matplotlib`](https://matplotlib.org)
- [`sklearn`](http://scikit-learn.org/stable/)
- [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) and the CPLEX Python API.
  Note: Do not install the Python CPLEX API by doing something like `pip install cplex`.
  This will install the trial version, which cannot handle large datasets.
  Instead, download CPLEX from the website, and run `python setup.py install` within the Python API folder.
  
## Subcommands

### Consensus Clustering
The input for consensus clustering is three files:
* A text file containing paths to distance matrices in `.tsv` format.
* A text file containing paths to clustering assignments in `.tsv` format.
* A text file containing the names of the clusterings which are 'finest'.

Each line should correspond to a specific data type, e.g. SNPs, MLSTs, or CNVs.
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

### Preprocessing
`MATCH_DIST` mode:
It may happen that for a particular pathogen, you do not have complete genotyping data for each sample, e.g. some samples may only have 1 or 2 out of the 3 genotyping datatypes, so that the distance matrices do not describe the same set of samples.
In this case, you cannot perform consensus clustering. 
The distance matrix matching subcommand takes as input multiple distance matrices and removes samples which don't appear in all of them so that the distance matrices describe the same set of samples.
We call this 'matching' the matrices.
The input for this mode is two files:
* A text file containing paths to the distance matrices.
* A text file containing output paths for the matched distance matrices. 
Each line should correspond to a specific data type, e.g. SNPs, MLSTs, or CNVs.
Paths to distance matrices should be prepended with the name of the datatype and an equal sign, i.e. `[name]=[path to file]`.
