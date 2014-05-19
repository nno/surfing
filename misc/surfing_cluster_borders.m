function b=surfing_cluster_borders(cl,nbrs)
% Finds the borders of clusters, useful for visualization
%
% B=SURFING_CLUSTER_BORDERS(CL,NBRS)
% 
% INPUT:
%   CL       Px1 cluster indices for P clusters; CL{K} are indices of
%            cluster K. Typically this is the output of SURFING_CLUSTERIZE
%   NBRS     KxN neigbours for K nodes, each with at most P neighbours
%            (values of zero denote no neighbour)
% OUTPUT:
%   B        Kx1 indices, where B(I)=K if CK{K}==I and I not on the border
%                                   =-K if CK{K}==I and I is on the border
%  
% Thus, to find indices that are on the border of a cluster, use find(B<0).
%
% Example: suppose x is a Nx1 data vector and v and f the coordinates and 
% faces for a surface. To find and visualize clusters with x>2, do:
%   n2a=surfing_surfacearea(v,f);
%   nbrs=surfing_surface_nbrs(f);
%   cl=surfing_clusterize(x>2,nbrs)
%   b=surfing_cluster_borders(cl,nbrs);
%   y=x;
%   y(b<0)=-100;
%   
% and now the borders of clusters are outlines by nodes with value -100
%
% See also SURFING_SURFACE_NBRS, SURFING_CLUSTERIZE
% 
% NNO Apr 2011

ncl=numel(cl); % number of clusters
[nv,mxnb]=size(nbrs); % vertices and max number of neighbours

% augment neighbours and nodes, by adding one dummy node at the very
% beginning. This esnures that all values in augnbrs are positive and thus
% indexing goes without issues
augnbrs=[repmat(1,1,mxnb); nbrs+1]; 
augnode2cl=ones(nv+1,1); % 

% mask to indicate if a node is on the border
onborder=false(nv,1);

% find which cluster each node is contained in
for k=1:ncl
    clk=cl{k};
    augnode2cl(clk+1)=k+1;
end

% find borders
for k=1:ncl
    clk=cl{k};  % indices of nodes in this cluster
    nbrsk=augnbrs(clk+1,:); % neigbours of all nodes
    nbrskn=sum(nbrsk>1,2); % number of neighbours of each node

    % find which cluster each neighbour is part of
    clknbrs=reshape(augnode2cl(nbrsk,:),[],mxnb);
    
    % node is on border if not all neighbours are part of this cluster
    bordermask=nbrskn>sum(clknbrs==k+1,2);
    
    % set border membership and update mask
    onborder(clk(bordermask))=true;
end

b=augnode2cl(2:end)-1;    % set cluster membership
b(onborder)=-b(onborder); % negative values for those 
