# SCANCell
SCANCell is a new analysis pipeline that unveils the hidden association networks of clusters for high dimensional cytometry data analysis. 

# Data preparation
SCANCell takes the data in .RData format as input. The CyTOF data of each sample is saved as a separate .fcs or .csv file. The data is firstly transformed by an arc-hyperbolic sine function (arcsinh(•)) with a cofactor of 5. Then data of all examples are concatenated into an integrated data of .RData format.

# Requirements
- R

- MATLAB

- Cytoscape

# Using SCANCell 
SCANCell mainly contains two analysis parts: sample classification and direct association (DA) network construction.

## Sample Classification
We use phEMD, proposed by Bodenmiller, to classify samples into different subgroups. The installing and running of phEMD refer to https://github.com/KrishnaswamyLab/phemd. We made slightly changes for phEMD and output the following results we required in a self-defined folder (such as the folder named results_from_PhEMD in the provided folder SCANCell).

- classification results of all samples

- cell clustering results of each sample

## DA Network Construction
For each subgroup of samples, taking their cell-marker matrix data and the corresponding clustering labels of cells as input, executing MATLAB code construct_network.m to obtain the DA network expressed as the edges with weight. During the network construction, we use part mutual information (see the manuscript Part mutual information for quantifying direct associations in networks by Zhao,J. et al.) to measure the direct interactions between clusters. Then the network is visualized by Cytoscape software. Some auxiliary comments are marked in construct_network.m file.

 
