function [nidxs,vidxs]=surfing_randomsubspace_select(c1,c2,voldef,nodeidxs,roisize,itercount,linedef,progressstep)
% selects voxels and nodes randomly from the surface, for a fixed number
% of voxels
%
% [NIDXS,VIDXS]=SURFING_RANDOMSUBSPACE_SELECT(C1,C2,VOLDEF,NODEIDXS,ROISIZE,ITERCOUNT,LINEDEF]
%
% INPUTS:
%   C1,C2      Px3 coordinates of inner and outer surface (white/pial)
%   VOLDEF     Struct with fields .dim and .mat (a la SPM_VOL)
%   NODEIDXS   1xQ indices of nodes of surface to select from
%   ROISIZE    number of voxels to select in each iteration
%   ITERCOUNT  number of iterations (i.e. ROIs)
%   LINEDEF    Default is [5 0 1] for 5 points from C1 to C2
% OUTPUTS:
%   NIDXS      1xITERCOUNT cell with node indices of selected voxels on surface
%   VIDXS      1xITERCOUNT cell with linear voxel indices of selected voxels. 
%              On average ROISIZE voxels indice are in each element
%
% NNO Oct 2010

if nargin<8 || isempty(progressstep), progressstep=100; end
if nargin<7 || isempty(linedef), linedef=[5 0 1]; end

% Find coordinates of points on or between the two surfaces.
% Here only for NODEIDXS (instead of whole surface)
allcoords=surfing_nodeidxs2coords(c1,c2,nodeidxs,linedef);

% Find the voxel indices corresponding to the coordinates of points above
% ALLLINVOXIDXS(I,K) contains the linear index for the voxel
% associated with node I for the K-th step. 
alllinvoxidxs=surfing_coords2linvoxelidxs(allcoords,voldef);

% For each row seperately, duplicates are replaced by zeros. 
unqlinvoxidxs=surfing_uniqueidxsperrow(alllinvoxidxs); 
clear alllinvoxidxs;

% construct intermediate surface that is the average of the two surfaces
intermediatecoords=squeeze(surfing_nodeidxs2coords(c1,c2,nodeidxs,[1,0.5,0.5]));

% find linear indices of voxels containing the center nodes. Values of NaN
% indicate outside the volume, and those are not used as center.
intermediatelinvoxidxs=surfing_coords2linvoxelidxs(intermediatecoords,voldef); 
outsidemsk=isnan(intermediatelinvoxidxs);
if sum(outsidemsk)>0 && progressstep    
    warning('surfing_voxelselection:note','found %d / %d center nodes outside the volume, these will be ignored.', sum(outsidemsk), ncenters);
end


nodecount=numel(nodeidxs);

% allocate space for output
nidxs=cell(itercount,1);
vidxs=cell(itercount,1);

% keep track of number of voxels and nodes
nsum=0;
vsum=0;

tic();

allunqcount=numel(unique([unqlinvoxidxs(:); 0]))-1;
if allunqcount<roisize
    warning('Cannot select %d indices: only %d voxels covered by %d nodes',roisize,allunqcount,numel(nodeidxs));
    nidxs=NaN;
    vidxs=NaN;
    return
end

for k=1:itercount
    rp=randperm(nodecount);
    rpvoxidxs=unqlinvoxidxs(rp,:);


    % indices of selected voxels, and number of nodes
    [sidxs,scount]=surfing_selectkfirstidxs(roisize,rpvoxidxs);
    
    nsum=nsum+scount;
    vsum=vsum+numel(sidxs);
    
    vidxs{k}=int32(sidxs);
    nidxs{k}=int32(rp(1:scount));
    
    if progressstep
        if k==1 || mod(k,abs(progressstep))==0 || k==itercount; % show progress in beginning, every PROGRESSSTEP nodes, and at the end
            if progressstep<0, clc(); end
            tc=toc();
            eta=(itercount-k)/k*tc;
            rtxt=sprintf('%.1f vx, %.1f nd',vsum/k,nsum/k);
            fprintf('Completed %d / %d nodes (%.1f%%), %s, took %d sec, ETA %d sec\n',k,itercount,k/itercount*100,rtxt,round(tc),round(eta));
        end
    end
end
    