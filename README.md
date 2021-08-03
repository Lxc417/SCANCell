# SCANCell
SCANCell is a new analysis pipeline that unveils the hidden association networks of clusters for high dimensional cytometry data analysis. SCANCell mainly contains two analysis parts: sample classification and direct association (DA) network construction.

# Data preparation
SCANCell takes the data in .RData format as input. The CyTOF data of each sample is saved as a separate .fcs or .csv file. The data is firstly transformed by an arc-hyperbolic sine function (arcsinh(•)) with a cofactor of 5. Then data of all examples are concatenated into an integrated data of .RData format.

# Requirements
R

MATLAB

Cytoscape

# Sample Classification
We use phEMD, proposed by Bodenmiller, to classify samples into different subgroups. The installing and running of phEMD refer to https://github.com/KrishnaswamyLab/phemd. We made slightly changes for phEMD and output the results we required in a self-defined folder.

classification results of all samples

cell clustering results of each sample

# DA Network Construction
For each subgroup of samples, taking their cell-marker matrix data and the corresponding clustering labels of cells as input, executing MATLAB code construct_network.m to obtain the DA network expressed as the edges with weight. Then the network is visualized by Cytoscape software. Some auxiliary comments are marked in construct_network.m file.

 
