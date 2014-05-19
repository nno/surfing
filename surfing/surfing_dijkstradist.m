function final_dist=surfing_dijkstradist(coords,faces,center_idx,max_distance,nbrs)
% Computes dijkstra distance on a surface mesh
%
% final_dist=surfing_dijkstradist(coords,faces,center_idx,max_distance,nbrs)
% 
% INPUTS:
%   coords:         3xP coordinates for P vertices
%   faces:          3xQ vertex indices for Q faces
%   center_idx:     vertex index from which the distances are computed
%   max_distance:   maximum distance to return (optional). Default is
%                   infinity.
%   nbrs:           PxM neighbor matrix, if each vertex has at most M
%                   neighbors. If omitted is computed based on faces.
%
% OUTPUT:
%   final_dist      Px1 distance to the vertex with index center_idx.
%                   Values of infinity are assigned if a node could not be
%                   reached or was further away than max_distance
%
% NNO Feb 2013

if nargin<4 || isempty(max_distance), max_distance=Inf; end

[three,nv]=size(coords);
[three_,nf]=size(faces);

if three~=3 || three_~=3, error('require 3x{V,F} faces and vertices'); end

% find the neighbors

if nargin<5 || isempty(nbrs)
    nbrs=surfing_surface_nbrs(faces');
end

% final distances
final_dist=Inf(nv,1);

% mask of candidates
candidates_mask=false(1,nv);
candidates_mask(center_idx)=true;

% tentative distances
tent_dist=Inf(1,nv);
tent_dist(center_idx)=0;

while sum(candidates_mask)>0
    candidates_idxs=find(candidates_mask);
    
    % find the candidate with the lowest tentative distance
    [d,i]=min(tent_dist(candidates_idxs));
    ifull=candidates_idxs(i);
    candidates_mask(ifull)=false; % not a candidate anymore
    
    if ~isinf(final_dist(ifull))
        continue % already has a final distance
    end
    
    % determine the neighbors
    nbr=nbrs(ifull,:);
    nbr=nbr(nbr>0);
    
    % take distance of candidate plus distance to each neighbor
    nbr_dist=d+surfing_eucldist(coords(:,ifull),coords(:,nbr)); 
    
    % see for which neighbors to update the tentative distance
    nbr_mask=nbr_dist<=max_distance & (isinf(tent_dist(nbr)) | nbr_dist<tent_dist(nbr));
    update_nbr_mask=nbr(nbr_mask);
    
    % update tentative distance
    tent_dist(update_nbr_mask)=nbr_dist(nbr_mask);
    candidates_mask(update_nbr_mask)=true;
    
    % set the final distance
    final_dist(ifull)=tent_dist(ifull);
end
    
    
    
    
    
    