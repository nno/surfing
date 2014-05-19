function [opts,b2s]=surfing_afni_opts_subvol(I,mn,mx)
% Options for loading minimal sub-volume with AFNI's BrikLoad
%
% [OPTS,ALL2FEW]=SURFING_AFNI_OPTS_SUBVOL(I,MN,MX)
% INPUTS:
%   I        AFNI Info struct (from BrikLoad), or 1x3 volume dimensions
%   MN,MX    List of minimum and maximum sub indices in x,y,z axes,
%            typically as returned by SURFING_VOXELSELECTION
% OUTPUTS
%   OPTS     options for AFNI's BrikLoad that load the minimal volume
%   ALL2FEW  Qx1 mapping of linear indices from 'entire' to 'minimal' volume
%
% Purpose: load minimal box of data necessary, to prevent out of memory
% issues
%
% Update (Jan 2010): With the new BrikLoad option '.Voxels', the present
% function is more or less obsolute. Use SURFING_REDUCEMAPPING instead,
% which allows for indicating for which voxels data is loaded.
%        
% NNO Oct 2010

if ischar(I)
    [err,I]=BrikInfo(brikfn);
end

if isstruct(I) && isfield(I,'DATASET_DIMENSIONS')
    dim=I.DATASET_DIMENSIONS(1:3);
elseif isstruct(I) && isfield(I,'dim')
    dim=I.dim;
elseif isnumeric(I)
    dim=I;
else
    error('Unrecognized struct I');
end

szmn=size(mn);
szmx=size(mx);
if ~isequal(szmn,szmx) || szmn(2) ~= 3
    error('mn and mx should be Px3');
end

mns=min(mn);
mxs=max(mx);

% range of voxel along x,y,z axes
xx=mns(1):mxs(1);
yy=mns(2):mxs(2);
zz=mns(3):mxs(3);

nx=numel(xx);
ny=numel(yy);
nz=numel(zz);
nxyz=nx*ny*nz;

% list all sub indices of voxels to include
subidxs=zeros(nxyz,3);
for i=1:nx
    for j=1:ny
        for k=1:nz
            subidxs((k-1)*nx*ny+(j-1)*nx+i,:)=[xx(i) yy(j) zz(k)];
        end
    end
end

linidxs=surfing_subs2inds(dim,subidxs);
b2s=zeros(prod(dim),1);
b2s(linidxs)=1:nxyz;

opts=struct();
opts.Format='vector';
opts.OutPrecision='*';
opts.Scale=0;
opts.PixX=xx;
opts.PixY=yy;
opts.Slices=zz;
