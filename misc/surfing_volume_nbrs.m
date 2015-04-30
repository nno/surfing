function nbrs=surfing_volume_nbrs(voldef,nn)
% Linear indices of neighbouring voxels in volume
%
% NBRS=SURFING_VOLUME_NBRS(VOLDEF,NN)
% INPUTS:
%   VOLDEF: spm volume definition from spm_vol, or AFNI filename or header
%           information, or 1x3 vector with number of voxels in x,y,z
%   NN      nearest neighbour level: any subset of [1,2,3] for touching
%           faces, edges, and/or corners. Default is [1,2,3]
% OUTPUT
%   NBRS:   PxQ neighbours, for P voxels and maximum of Q neighbours each.
%           Contains NaNs to indicate no neighbour.
%
% NNO Nov 2010


if nargin<2
    nn=[1 2 3];
elseif setdiff(nn,1:3)
    error('Illegal NN definition: should  be a subset of [1 2 3]')
end

if (ischar(voldef) && regexp(voldef,'+orig|+tlrc|+acpc')) ...
        || isstruct(voldef) && isfield(voldef,'DATASET_DIMENSIONS')
    % AFNI file?
    voldef=surfing_afni2spmvol(voldef);
end

if isstruct(voldef) && isfield(voldef,'dim')
    dim=voldef.dim;
elseif isnumeric(voldef) && numel(voldef)==3
    dim=voldef;
else
    error('Unrecognized voldef');
end

ndim=numel(dim);

offsets=surfing_inds2subs(repmat(3,1,ndim),(1:(3^ndim))')-2; % 27x3 matrix with offsets relative to center voxel
offsetsnn=sum(offsets.^2,2); % NN level for offsets
nnnbrs=histc(offsetsnn,1:ndim); % should be [6, 12, 8]'

nnmsk=false(1,ndim); % NN 1,2,3?
nbrcount=0; % max number of neighbours
for k=1:ndim
    addk=sum(nn==k)>0;
    nnmsk(k)=addk;
    nbrcount=nbrcount+addk*nnnbrs(k);
end

voxcount=prod(dim); % number of voxels in volume
ijkpos=surfing_inds2subs(dim,(1:voxcount)'); % ijk indices of voxels in volume

nbrs=NaN(voxcount,nbrcount); % allocate space for output
freecol=0; % first free column in output
for k=1:ndim
    if ~nnmsk(k)
        continue;
    end

    offsetidxs=find((k==offsetsnn));
    for j=1:numel(offsetidxs)
        offsetsj=offsets(offsetidxs(j),:);
        ijkrel=bsxfun(@plus,offsetsj,ijkpos);
        linrel=surfing_subs2inds(dim,ijkrel);
        nbrs(:,freecol+j)=linrel;
    end

    freecol=freecol+nnnbrs(k);
end



