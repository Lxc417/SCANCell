clear;
labels=[];
Data=[];
AllFiles=dir('data\*.fcs');
%% reading data
for i=1:numel(AllFiles)
    ll=readtable(['results_from_PhEMD\labels_',AllFiles(i).name,'.csv']); %read the cell clustering label of each sample in a subgroup
    l_l=ll.x;
    [fcsdat, fcshdr]=fca_readfcs(['data\',AllFiles(i).name]); %read the cell-marker data of each sample in a subgroup
    labels=[labels;l_l]; %merge all the clustering labels of a subgroup
    Data=[Data;fcsdat]; %merge all the cell-marker data of a subgroup   
end

markers={'CD5','CD8','CD1d','CD185','CD279','CD56','CD183','FOXP3','RORyt','CD25','CD3','CD27','CD4','CD196'};

%% calculating the median cluster-marker matrix for clusters
[n,m]=size(Data);
C_size = length(unique(labels)); %the number of clusters
ClassIndex=cell(1,C_size);
X=cell(1,C_size);
X_median=zeros(C_size,m);
for i=1:C_size
    ClassIndex{i}=find(labels==i); %the cell index of each cluster
    clusterNum(i)=length(ClassIndex{i}); %the number of cells of each cluster
    X{i}=Data(ClassIndex{i},:); %extract the cell-marker data for each cluster
    X_median(i,:)=median(X{i},1); %calculate the median expression for each cluster
end

for i=1:C_size
    clusterName(i)={['C-',num2str(i)]};
end
cell_freq=clusterNum'/sum(clusterNum); %the cell abundance of each cluster

bar(cell_freq);
title('UN\_RP')
xlabel('cell clusters','FontSize',12)
ylabel('cell frequency','FontSize',12)
axis([0 C_size+1 0 max(cell_freq)])
set(gca,'FontSize',12);
set(gca,'xtick',1:C_size,'xticklabel',clusterName) 
xtickangle(45)

figure;
imagesc(X_median)
colorbar
colormap(jet)
title('cell clusters - markers:UN\_RP')
xlabel('markers','FontSize',12)
ylabel('cell clusters','FontSize',12)
axis([0.5 m+0.5 0.5 C_size+0.5])
set(gca,'FontSize',12);
set(gca,'ytick',1:C_size,'yticklabel',clusterName)
set(gca,'xtick',1:m,'xticklabel',markers)
xtickangle(45)

%% Calculating PMI to measure the direct association between two clusters and saving network data in .csv files
[G,Gval,l]=pca_pmi(X_median,2,10); %call pca_pmi function to calculate the PMI values of any paired clusters.In PMI algorithm, a threshold of reducing edges is taken as 2 and the processing order is set to 10. When the processing order is larger than 10, the algorithm will stop reducing the edge and return the network.
    
weight=Gval+Gval';
A=G+G';
d=sum(A)';
[row,col]=find(G~=0);
    
val=nonzeros(Gval);
val=val(val > 2);
NE=[row col val];

var={'node1','node2','weight'};
T=array2table(NE,'VariableNames',var);
writetable(T,'net.csv') %save network data represented by edges and corresponding weights
dlmwrite('weight.csv',Gval) %save the weighted adjacency matrix


figure;
imagesc(weight)
colorbar
colormap(jet)
title('interactive weight:UN\_RP')
xlabel('cell clusters','FontSize',12)
ylabel('cell clusters','FontSize',12)
axis([0.5 C_size+0.5 0.5 C_size+0.5])
set(gca,'FontSize',12);
set(gca,'xtick',1:C_size,'xticklabel',clusterName)
set(gca,'ytick',1:C_size,'yticklabel',clusterName)
xtickangle(45)




