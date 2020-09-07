clear;
labels=[];
Data=[];
AllFiles=dir('data\*.fcs');
%% reading data
for i=1:numel(AllFiles)
    ll=readtable(['results_from_PhEMD\labels_',AllFiles(i).name,'.csv']);
    l_l=ll.x;
    [fcsdat, fcshdr]=fca_readfcs(['data\',AllFiles(i).name]);
    labels=[labels;l_l];
    Data=[Data;fcsdat];    
end

markers={'CD5','CD8','CD1d','CD185','CD279','CD56','CD183','FOXP3','RORyt','CD25','CD3','CD27','CD4','CD196'};

%% calculating the median cluster-marker matrix for clusters
[n,m]=size(Data);
C_size = length(unique(labels));
ClassIndex=cell(1,C_size);
X=cell(1,C_size);
X_median=zeros(C_size,m);
for i=1:C_size
    ClassIndex{i}=find(labels==i);
    clusterNum(i)=length(ClassIndex{i});
    X{i}=Data(ClassIndex{i},:);
    X_median(i,:)=median(X{i},1);
end

for i=1:C_size
    clusterName(i)={['C-',num2str(i)]};
end
cell_freq=clusterNum'/sum(clusterNum);

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

%% PMI calculation and saving network data
[G,Gval,l]=pca_pmi(X_median,2,10);
    
weight=Gval+Gval';
A=G+G';
d=sum(A)';
[row,col]=find(G~=0);
    
val=nonzeros(Gval);
val=val(val > 2);
NE=[row col val];

var={'node1','node2','weight'};
T=array2table(NE,'VariableNames',var);
writetable(T,'net.csv')
dlmwrite('weight.csv',Gval)


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




