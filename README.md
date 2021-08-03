# SCANCell
SCANCell is a new analysis pipeline that unveils the hidden association networks of clusters for high dimensional cytometry data analysis. 
#Data preparation
SCANCell takes the data in .RData format as input. The CyTOF data of each sample is saved as a separate .fcs or .csv file. The data is firstly transformed by an arc-hyperbolic sine function (arcsinh(•)) with a cofactor of 5. Then data of all examples are concatenated into an integrated data of .RData format.

SCANCell first depicts the immune phenotype of  biological cohorts and classifies samples from the same experimental group  into different subgroups using PhEMD (https://github.com/KrishnaswamyLab/phemd). 
Then SCANCell  computes the PMI values between clusters and constructs DA networks. 
