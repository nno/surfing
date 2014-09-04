function [n2ns,radii]=surfing_nodeselection(v,f,circledef,dist_metric,progressstep)
% gives indices of neighboring nodes for each node
%
% [n2ns,ds]=surfing_nodeselection(v,f,circledef[,dist_metric])
% 
% Inputs:
%   v               3xP node coordinates for P nodes
%   f               3xQ face indices for Q faces
%   circledef       searchlight radius (in mm), or NaN to return direct
%                   neighbors only
%   dist_metric     'euclidean', 'dijkstra', or 'geodesic' (default)
%
% Outputs:
%   n2ns            Px1 cell, with n2ns{k} containing the indices of nodes
%                   within distance circledef from node k according to the
%                   specified distance metric
%   ds              Px1 radius of each searchlight
%
% NNO June 2014

if size(v,1)~=3 || size(f,1) ~= 3
    error('vertices and faces must have 3 elements in first dimension');
end

if nargin<4, dist_metric='geodesic'; end
if nargin<5 || isempty(progressstep), progressstep=100; end

show_progress=~isempty(progressstep) && progressstep~=0;

nv=size(v,2);
n2ns=cell(nv,1); % mapping from center node to surrounding nodes

if isnan(circledef)
    [n2ns_matrix,all_radii]=surfing_surface_nbrs(f',v');
    for k=1:nv
        node_indices=n2ns_matrix(k,:);
        n2ns{k}=node_indices(node_indices>0);
    end
    radii=max(all_radii,[],2);
    return
end

% for speedup, precompute mapping 
n2f=surfing_invertmapping(f');


% space for output
radii=zeros(nv,1);
sizes=zeros(nv,1);

clock_start=clock();
prev_msg='';

for k=1:nv
    [idxs,distances]=surfing_circleROI(v,f,k,circledef,dist_metric,n2f);
    n2ns{k}=idxs;
    sizes(k)=numel(idxs);
    if isempty(distances)
        distances=Inf;
    end
    radii(k)=max(distances);
    
    if show_progress && (k<10 || mod(k,progressstep)==0 || k==nv)
        msg=sprintf('%.1f nodes per center', mean(sizes(1:k)));
        prev_msg=surfing_timeremaining(clock_start,k/nv,msg,prev_msg);
    end
end
    