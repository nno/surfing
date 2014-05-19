function [b2s,s2b,x2ys]=surfing_reducemapping(x2y,ymask)
% simplifies a cell mapping
%
% [BIG2SMALL,SMALL2BIG,X2YSMALL]=SURFING_REDUCEMAPPING(X2Y[,MASK])
% INPUT:
%   X2Y          Nx1 cell, each with a vector with positive integers
%   MASK         optional mask indicating which integers in X2Y should be
%                excluded from the mapping in the output
% OUTPUT:
%   BIG2SMALL    Mx1 vector that maps all integers in X2Y to the indices
%                of the unique integers in X2Y. M is the maximum value
%                across all vectors in X2Y. Non-mapped integers have a
%                value of zero.
%   SMALL2BIG    Ux1 The unique integers in X2Y in ascending order
%   X2YSMALL     Nx1 cell, each with a vector with positive integers in the
%                range 1:U.
%
% It holds that for all i: SMALL2BIG(BIG2SMALL(X2Y{i})) == X2Y{i}
%                 ... and: BIG2SMALL(X2Y{i})            == X2YSMALL{i}
% 
% Example: suppose X2Y is a mapping from center nodes to surrounding
% voxel indices, referring to a matrix BETA where each row corresponds to
% values associated with a single voxel. Many voxels may not be associated 
% with any voxel, for example if they are outside the brain.
% Assuming IDX is the index of a center node, then BETA(X2Y{IDX}, :) are 
% the voxel values of the surrounding nodes.
%
% To reduce memory, use [BIG2SMALL,SMALL2BIG]=SURFING_REDUCEMAPPING(X2Y)
% followed by BETASMALL=BETA(SMALL2BIG,:).  (BETA can be cleared from memory) 
%
% The voxel values associated with IDX can now be accessed through
% BETASMALL(BIG2SMALL(X2Y{IDX}),:)
%
% -----------------------------------------------------------------------
%
% For AFNI users: as of Mar 2011, the AFNI matlab function BrikLoad supports
% an argument opt.Voxels, so that for information mapping purposes loading 
% functional data can be restricted to voxels that are in one or more
% searchlights. In the following example, n2v is the output from 
% surfing_voxelselection, epibrikfn is the file with functional data, and 
% fnout is the surface .niml.dset file name to which output is written
%   
%   [b2s,unq,n2vs]=surfing_reducemapping(n2v);
%   opt=struct(); opt.Voxels=unq; opt.Format='vector';
%   [err,V,I]=BrikLoad(epibrikfn,opt); % load functional data
%
%   n=numel(n2vs); % number of searchlights
%   r=zeros(n,1); % output result
%   for k=1:n
%       Vk=V(n2vs{k},:); % data in k-th searchlight
%       % ... do something with Vk ...
%       r(k)=...
%   end
%
%   S=struct(); S.data=r; S.node_indices=unq-1; % output
%   afni_niml_writesimple(S,fnout);
%
% NNO Dec 2010

isvec=isnumeric(x2y); % allow for just a single vector as input
if isvec
    x2yvec=x2y(:)';
    x2y=cell(1);
    x2y{1}=x2yvec;
end

if ~iscell(x2y), error('x2y should be a cell'); end

n=numel(x2y);

% find maximum value across all cell vectors
maxy=0;
for k=1:n
    x2yk=x2y{k};
    if isempty(x2yk)
        continue;
    elseif ~all(isfinite(x2yk)) || ~isequal(round(x2yk),x2yk) || min(x2yk)<1
        error('vector at x2y{%d} does not contain only positive integers',k);
    end
    
    m=max(x2y{k});
    if m>maxy
        maxy=m;
    end
end

% define a mask of values found in any of x2y
msk=false(maxy,1); % hopefull maxy is not too big... otherwise we run out of memory
for k=1:n
    msk(x2y{k})=true; % positive integers assumed in x2y{k}
end

% if mask is provided, apply it now
if nargin>=2
    nymask=numel(ymask);
    
    showwarning=true;
    if isnumeric(ymask) % if numeric, convert to boolean; warning off
        ymask=find(ymask);
        showwarning=false;
    end
    
    if islogical(ymask)
        if showwarning && nymask~=maxy
            warning('reducemapping::number of elements in mask for y (%d) does not match maximum of x2y (%d)',nymask,maxy);
        end
        maxn=min(numel(ymask),maxy);
        
        msk(1:maxn)=msk(1:maxn) & ymask(1:maxn); % intersect 
    end
end
    
s2b=find(msk);
m=max(s2b);
b2s=zeros(1,m);
b2s(s2b)=1:numel(s2b);

if nargout>=3
    x2ys=cell(size(x2y));
    for k=1:n
        x2ys{k}=b2s(x2y{k});
    end
    
    if isvec
        x2ys=x2ys{1};
    end
end