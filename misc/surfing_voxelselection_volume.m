function [v2vs,mn,mx,vORr]=surfing_voxelselection_volume(spheredef,voldef,centeridxs,progressstep)
% volume based voxel selection for searchlight
%
% [V2VS,MN,MX,VORR]=SURING_VOXELSELECTION_VOLUME(SPHEREDEF,VOLDEF[,CENTERIDXS[,PROGRESSSTEP]])
%
% INPUTS:
%   SPHEREDEF:   definition of searchlight sphere, either one of:
%                [R]   scalar R means fixed radius R (in mm)
%                [R C] use maximum radius R to find approximately C voxels
%   VOLDEF:      volume definition struct with fields:
%                .dim       1x3 vector with volume dimensions (in voxels)
%                .voxsize   1x3 vector with voxels size (in mm)
%                .mask      vector or array with PROD(VOLDEF.voxsize)
%                           elements (logical or numerical) with brain mask
%                For AFNI users, use SURFING_AFNI2SPMVOL to get a voldef
%   CENTERIDXS:  1xP vector with linear indices of searchlight center voxels
%                all values should be in range 1:PROD(VOLDEF.dim)
%                If omitted, default is all center indices
%   PROGRESSTEP  show progress every PROGRESSTEP voxels 
%                (default: 1000, use 0 for mute)
% OUTPUTS:
%   V2VS         1xP cell, with V2VS{k} containing the linear indices of the
%                voxels surrounding the voxel with index CENTERIDXS(k)
%   MN,MX        Px3 vectors with minimum and maximum voxel sub indices
%   VORR         1xP vector with searchlight radii or number of voxels
%
% This function is not a nice as the surface equivalent:
% (1) requires non-oblique functional volumes
% (2) no automatic increase of sphere size if it is too small
%
% NNO Sep 2010

% process spheredef
sphereradius=spheredef(1);
fixedradius=numel(spheredef)==1;
if ~fixedradius
    targetvoxelcount=spheredef(2);
end

% process voldef
dim=voldef.dim;
voxcount=prod(dim);

if isfield(voldef,'mat') % support affine matrix (non oblique)
    mat=voldef.mat;
    mat(4,4)=1;
    voxsize=[mat(1,1) mat(2,2) mat(3,3)];
    mat=mat-diag([voxsize 1]);
    if ~isequal(mat(1:3,1:3),zeros(3))
        % maybe later we support oblique matrices, but not for now
        kk=mat(1:3,1:3);
        ukk=unique(kk(kk~=0));
        ukk
        if numel(ukk)==1
            voxsize=repmat(ukk,1,3);
        else
            voldef.mat
            error('matrix oblique, not supported yet');
        end
    end
elseif isfield(voldef,'voxsize') % alternative, use voxsize
    voxsize=voldef.voxsize;
    mat=diag([voldef.voxsize 1]);
else
    error('expected voldef with field .mat or .voxsize');
end

if isfield(voldef,'mask') % see if there is a mask
    mask=voldef.mask(:);
    if numel(mask) ~= voxcount % check size
        error('size of voldef.mask and prod(voldef.dim) do not match');
    end
    if isnumeric(mask)
        mask=logical(mask);
    end
else
    mask=true(voxcount,1); % use all voxels
end

if nargin<3 || isempty(centeridxs)
    centeridxs=1:voxcount; % use all voxels as center
else
    allcenteridxscount=numel(centeridxs);
    selectidxscount=numel(intersect(centeridxs,find(mask)));
    warning('surfing:voxelselection_volume','%d / %d (%.1f %%) center voxels outside the mask!',allcenteridxscount-selectidxscount,allcenteridxscount,(allcenteridxscount-selectidxscount)*100/allcenteridxscount);
end

if nargin<4 || isempty(progresstep)
    progressstep=1000; % default: show progress every 1000 voxels
end

centercount=numel(centeridxs); % number of center voxels

centerijks=phoebe_inds2subs(dim,centeridxs(:)); % sub indices of center voxels
voxorder=randperm(centercount); % permute visiting order to get good estimate of number of voxels selected

% offsets for spherical searchlight mask
[maskoffsets,allcenterdistances,spherevolmask]=phoebe_spherical_mask(sphereradius,abs(voxsize));
centerdistances=allcenterdistances(spherevolmask(:));

lin2sub=surfing_inds2subs(voldef.dim,1:voxcount);


% some variables during voxelselection
v2vs=cell(1,centercount); % allocate space for output
vORr=zeros(centercount,1);
rs=zeros(centercount,1);
ascenter=false(centercount,1); % which voxels were used as center
aroundcountsum=0; % how many voxels selected total
mn=NaN(voxcount,3); % minimum and maximum voxel indices
mx=NaN(voxcount,3);



tic(); % start the clock
for k=1:centercount % loop of center voxels
    voxk=voxorder(k);
    centeridx=centeridxs(voxk);
    ignorevoxel=~mask(centeridx); % if outside the mask, ignore this as the center
    
    aroundcount=0;
    if ~ignorevoxel
        centerijk=centerijks(voxk,:); % sub index for this center voxels
        aroundijks=bsxfun(@plus,centerijk,maskoffsets); % add offsets to get sub indices of surrounding voxels
        aroundidxs=phoebe_subs2inds(dim,aroundijks); % convert to linear indices
        
        keepmask=~isnan(aroundidxs); % see which voxels are inside the volume
        
        aroundidxskeep=aroundidxs(keepmask); % indices around center node inside volume
        aroundidxsmask=mask(aroundidxskeep); % get mask values for the voxels to keep 
        
        keepmask(keepmask)=aroundidxsmask;  
        
        if fixedradius
            % easy!
            radiusk=sphereradius;
        else
            % more tricky: find the distances of the voxels we have kept
            % and select the nearest one
            arounddistances=centerdistances(keepmask);
            [foo,sidxs]=sort(arounddistances);

            lastidx=min(sum(keepmask),targetvoxelcount);
            tofarmask=true(sum(keepmask),1);
            tofarmask(1:lastidx)=false;
            
            keepmask(keepmask)=~tofarmask;
            radiusk=arounddistances(sidxs(lastidx));
        end
        aroundidxs=aroundidxs(keepmask); % also get rid of voxels outside the mask
        
        aroundcount=numel(aroundidxs); % number of voxels in searchlight
        if aroundcount==0
            ignorevoxel=true;
        elseif ~fixedradius && aroundcount~=targetvoxelcount
            warning('Voxel %d: selected only %d < %d voxels', voxk, aroundcount, targetvoxelcount);
        end
    end
    
    aroundcountsum=aroundcountsum+aroundcount;
    ignorevoxel=ignorevoxel || aroundcount==0;
    
    if ~ignorevoxel
        v2vs{voxk}=uint32(aroundidxs)'; % save linear indices of voxels in searchlight
        mn(voxk,:)=min(lin2sub(aroundidxs));
        mx(voxk,:)=max(lin2sub(aroundidxs));
        
        if fixedradius
            vORr(voxk)=aroundcount;
        else
            vORr(voxk)=radiusk;
        end
        ascenter(voxk)=true; % mark as center
        rs(voxk)=radiusk;
    end
    
    if progressstep % show progress?
        if k==1 || mod(k,abs(progressstep))==0 || k==centercount; % show progress in beginning, every PROGRESSSTEP nodes, and at the end

            if progressstep<0, clc(); end
            tc=toc();
            eta=(centercount-k)/k*tc;
            
            ascentercount=sum(ascenter); 
            meanradius=mean(rs(ascenter));
            rtxt=sprintf('r=%.2f, %.1f vox',meanradius,aroundcountsum/ascentercount); % radius and voxel count text
            
            fprintf('Completed %d / %d center voxels (%.1f%%), %s, took %d sec, ETA %d sec\n',k,centercount,k/centercount*100,rtxt,round(tc),round(eta));
        end
    end
    
end
    
