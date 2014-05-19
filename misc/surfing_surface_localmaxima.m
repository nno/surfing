function [is,mx]=surfing_surface_localmaxima(v,f,data,mindist,thr)
% Finds nodes that are a local maximum on the surface
%
% IDXS,MX=SURFING_LOCALMAXIMA(V,F,DATA,MINDIST,THR)
% INPUTS:
%   V        3xP surface coordinates for P nodes
%   F        3xQ node indices for Q triangles 
%   DATA     1xP data on which maxima are to be found
%   MINDIST  Minimum distance between local maxima
%   THR      Treshold to which data is submitted to before maxima are found
% OUTPUTS:
%   IDXS     Indices of nodes with local maxima n descending order
%   MX       Local maxima values
%
% Distances are measured with geodesic metric. Nodes returned in IDXS have 
% at least value THR and are at least MINDIST apart. 
%
% In other words: for each node with index IDX in IDXS, all nodes within 
% distance MINDIST from V(:,MX(IDX), all other nodes have a value less than
% or equal to MX(IDX)
%
% If THR is negative, then this function returns local minima of nodes that
% are less than THR
%
% NNO Sep 2011
% Earlier version posted online (23/4/2011):
% http://afni.nimh.nih.gov/afni/community/board/read.php?f=1&i=37971&t=37941&v=f
if thr<0
    thr=-thr;
    data=-data;
end

nbrs=surfing_surface_nbrs(f); % find neighbours
[nv,maxnb]=size(nbrs); % number of vertices, and max number of neighbours

if nv~=numel(data), error('Data should have %d values', nv); end

localmax=false(nv,1); % whether each node is a local maximum

for k=1:nv
    if data(k)<thr % skip nodes not exceeding threshold
        continue;
    end
    ismax=true;
% see if all neighbours are smaller than the value in data(k)
    for j=1:maxnb
        if nbrs(k,j)==0
            continue
        elseif data(k) <= data(nbrs(k,j))
            ismax=false;
            break
        end
    end
    if ismax
        localmax(k)=true;
    end
end

localmaxidxs=find(localmax);
nlocalmax=numel(localmaxidxs);

fprintf('Initially %d / %d nodes have a local maximum\n', nlocalmax, nv);

% now remove some nodes if they do not meet min_dist criteria
[mx,is]=sort(data(localmax),'descend'); % sort local maxima for min_dist
n2f=surfing_nodeidxs2faceidxs(f'); % find faces of each node, for faster operation

% loop over nodes (from large to small), remove all those that are within
% mindist geodesic distance
tic();
for k=1:nlocalmax
    maxidx=localmaxidxs(is(k));
    if ~localmax(maxidx)
        continue; 
    end % skip, not a local maximum anymore
    aroundidxs=surfing_circleROI(v',f',maxidx,mindist,n2f);
    localmax(aroundidxs)=false; % remove other nearby local maxima
    localmax(maxidx)=true; % but retain the current node
    surfing_timeremaining(k/nlocalmax);
end

localmaxidxs=find(localmax);
nlocalmax=numel(localmaxidxs);

fprintf('After applying mindist=%.1fmm, %d / %d nodes have a local maximum\n', mindist, nlocalmax, nv);

[mx,iso]=sort(data(localmax),'descend'); % sort again

is=localmaxidxs(iso);
