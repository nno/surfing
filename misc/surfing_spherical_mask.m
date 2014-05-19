function [ijkoffsets,distbox,m3]=surfing_spherical_mask(radius,voxelsize)
% constructs a sphere-shaped mask for use in a 3D volume (e.g. searchlight)
%
% [OFFSETS,BOX,M]=SURFING_SPHERICAL_MASK(R,VOXSIZE) returns 
% a set of linear offset values and a box with computed distances from the 
% center, using a mask radius R and a volume with voxelsize VOXSIZE.
% R should be a single scalar, and VOXSIZE a 1 x 3 vector.
% A binary mask M is also returned
%
% [...]=SURFING_SPHERICAL_MASK(R,INFO) uses a INFO struct 
% as obtained by AFNIs BrikInfo, and reads voxel size from its
% .DELTA 
%
%
% NNO Nov 2009

if nargin==2 && isstruct(voxelsize) && isfield(voxelsize,'DELTA')
    voxelsize=abs(voxelsize.DELTA);
end

[one,dimcount]=size(voxelsize); % in an earlier implementation, dimcount could be any number, but support for only dimcount=3 seems sensible
if one ~= 1 || dimcount ~= 3
    error('voxel size should be 1x3');
end

if numel(radius) ~= 1 || radius <= 0
    error('radius should be a positive scalar')
end

% required radius, in voxel sized units
halfsizes=floor(repmat(radius,1,dimcount)./voxelsize);

% box in which distances are stored
distbox=zeros(2*halfsizes+1);
boxsize=size(distbox);

% linear indices for each element in the box
boxcount=prod(boxsize);
linidxs=(1:boxcount)';

% convert to ijk coordinates 
ijkidxs=phoebe_inds2subs(boxsize,linidxs);

% center voxel
%centerijk=halfsizes+1;
centerlinidx=(boxcount+1)/2;

% position of each index, in world space
xyzcoords=ijkidxs .* repmat(voxelsize,boxcount,1);
distances=surfing_eucldist(xyzcoords(centerlinidx,:)',xyzcoords');

mask=distances<=radius;
distbox(:)=distances;
ijkoffsets_all=ijkidxs-repmat(ijkidxs(centerlinidx,:),boxcount,1);
ijkoffsets=ijkoffsets_all(mask,:);

m3=distbox<=radius;