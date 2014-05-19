function [nodeidxs, minds]= surfing_voxel2node(vxyz, nxyz)
% maps voxels coordinates to node coordinates 
%
% [NODEIDXS, MINDS]=SURFING_VOXEL2NODE(VXYZ,NXYZ) 
% INPUTS:
%   VXYZ:       3xP matrix for P voxel coordinates
%   NXYZ:       3xQ matrix for Q node coordinates
% OUTPUT
%   NODEIDXS:   node index for each voxel coordinate VXYZ
%   MINDS:      the distance from the voxel to the corresponding node
%
% TW July 2010

[three1,n1]=size(vxyz);
if three1 ~= 3, error('xs should be 3xN'); end

[three2,n2]=size(nxyz);
if three2 ~= 3, error('ys should be 3xN'); end


%----calc distances
ds=surfing_eucldist(vxyz,nxyz); % voxelxnodesCoord matrix
%----find the closest node
[minds, nodeidxs]= min(ds,[],2);
%----if the minimal distances are getting to big warn the user
if max(minds)>4
    warning('voxels have a distance that is bigger than 4mm!') 
end

