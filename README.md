# Node-Centric Walk Kernels
Source code of the paper [Weisfeiler and Leman Go Walking: Random Walk Kernels Revisited](https://doi.org/10.48550/arXiv.2205.10914), Nils M. Kriege, NeurIPS 2022.

## Usage
The graph kernels contained in this package can be computed via a command line interface. Run the shell script `kkernel` to see a list of all available kernels and parameters.

### Example
The following command computes the node-centric $\ell$-walk kernel for $\alpha$ in $\{0.1, 1.0\}$, $\beta$ in $\{0.0, 1.0\}$ and length $\ell$ in $\{0, 1, 2, 3, 4, 5\$} for the data set MUTAG:
```
./kkernel -d MUTAG ncw -a 0.1, 1.0 -b 0.0, 1.0 -l 0, 1, 2, 3, 4, 5
```
For each parameter combination the kernel matrix is computed and stored in the directory `gram` using the [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) file format.

## Building from source
Run `ant` to build `kgraph.jar` from source. 

## Data sets
The repository contains the data set MUTAG only. Further data sets are available from the [TUDataset](http://graphlearning.io) website. Please note that in the experimental comparison the edge labels, if present, were ignored. In order to reproduce the published results, please delete the files `DS_edge_labels.txt`, where `DS` is the name of the data set.

## Contact information
If you have any questions, please contact [Nils Kriege](https://dm.cs.univie.ac.at/team/person/111520/#info).
