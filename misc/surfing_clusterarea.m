function cl2area=surfing_clusterarea(clu,node2area)
% quick wrapper function to find areas of clusters. 
%
% CL2AREA=SURFING_CLUSTERAREA(CL,NODE2AREA)
%
% INPUTS:
%   CL         Nx1 cell for N clusters, where CL{K} are the indices of 
%              the nodes in the K-th cluster
%   NODE2AREA  Mx1 areas of each node
% OUTPUT:
%   CL2AREA    Nx1 areas for each cluster
%
% NNO Dec 2011


if ~iscell(clu), error('Clu should be a cell'); end
    
nclu=numel(clu);
cl2area=zeros(nclu,1);
for k=1:nclu
    cluidxs=clu{k};
    cl2area(k)=sum(node2area(cluidxs));
end
    