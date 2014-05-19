function V=surfing_afni2spmvol(I)
% returns an affine matrix for an AFNI brik similar to SPM
%
% V=SURFING_AFNI_AFFINEMATRIX(INFO)
% INPUT:
%   INFO:  Info brik as returned by BrikInfo (or a filename of an AFNI brik)
% OUTPUT:
%   V:      struct with the following fields:
%   V.mat:   4x4 affine transformation matrix, LPI convention (a la SPM)
%   V.dim: 1x3 vector with number of voxels in x,y,z
%
% This function requires the AFNI matlab toolbox
%
% NNO June 2010
%
% See also BRIKINFO

if isstruct(I) && isfield(I,'dim')
    V=I; %seems like an SPM struct; we are done
    return
end

if ischar(I)
    [err,I]=BrikInfo(I);
    if err
        error('Error loading %s', I);
    end
end

orient='LPI'; % always return LPI-based matrix

% origin and basis vectors in world space
k=[0 0 0;eye(3)];

[err,i]=AFNI_Index2XYZcontinuous(k,I,orient);

% basis vectors in voxel space
e1=i(2,:)-i(1,:);
e2=i(3,:)-i(1,:);
e3=i(4,:)-i(1,:);

% change from base0 (afni) to base1 (SPM/Matlab)
o=i(1,:)-(e1+e2+e3);

% create matrix
k=[e1;e2;e3;o]';

% set 4th row
k(4,:)=[0 0 0 1];

% set dimensions
V=struct();
V.mat=k;
V.dim=I.DATASET_DIMENSIONS(1:3);
V.fname=I.RootName;