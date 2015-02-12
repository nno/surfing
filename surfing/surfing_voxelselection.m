function [n2v,voxmin,voxmax,vORr]=surfing_voxelselection(c1,c2,f,circledef,voldef,nodeidxs,linedef,distancemetric,progressstep)
% Voxelselection function for a surface-based circular searchlight
% performes voxelselection using two surfaces and a circular ROI 
% 
% [N2V,VMIN,VMAX,DS]=SURFING_VOXELSELECTION(C1,C2,F,CIRCLEDEF,VOLDEF,[,NODEIDXS,LINEDEF,DIST,STEP])
% INPUT: 
%   C1          3xP coordinates for P vertices
%   C2          either:
%               - 3xP coordinates for P vertices; in this case C1 and C2
%                 are typically the pial and white surface from freesurfer
%               - 1x2 vector [offset_start, offset_stop] indicating the
%                 distances 'inwards' and 'outwards' along the normal 
%                 vectors of C1 along which a strip of brain is delineated;
%                 in this case C1 is usually a surface from brainvoyager or
%                 caret. For example, [-2 3] would select voxels along the
%                 normals of C1 between -2 units (typically mm) inwards and
%                 3 units outwards.
%   F:          3xQ vertex indices for Q faces (1-based indices)
%   CIRCLEDEF:  - scalar R, or [R 0]: for fixed searchlight radius R, 
%                                     variable number of voxels per node. 
%               - vector [R N]: fixed number of N voxels per node, 
%                               variable radius (initially radius R, use
%                               R=0 for auto)
%   VOLDEF:     Struct with fields .mat (4x4 voxel to coordinate transformation matrix 
%               from the images to be sampled (1-based)) and .dim (1x3
%               volume dimensions in voxels). Optionally a field .mask with
%               a voxel mask 
% OPTIONAL INPUTS:
%   NODEIDXS:   Sx1 vector for the center vertices for ROIS ([] for all, meaning S==P)
%   LINEDEF:    definition for lines spanning from surface 1 to surface 2 
%               (see SURFING_COORDS2VOXELIDXS)
%   DIST:       Distance metric 'geodesic' (default) or 'euclidian' or
%               'dijkstra'
%   STEP:       Progress reporting every STEP center vertices (0 for mute) 
% OUTPUT: 
%   N2V:        Sx1 cell so that N2V{K} is a 1xM_K vector of type uint32, containing 
%               the M_K voxel linear indices of the selected voxels around vertex K.
%               Use IND2SUB or SURFING_INDS2SUBS(DIM,N2V{K}) to get the subindices.
%   VMIN,VMAX:  Sx3 matrices with minimum and maximum voxel indices (i,j,k)
%   DS:         Sx1 vector with either number of voxels selected (for fixed radius),
%               or the radius for the searchlight (for fixed number of
%               voxels)
%
% - Distances are measured along the surface that is the vertex-wise average of 
%   the coordinates in C1 and C2. 
% - The function removes dublicates of voxel indices
% - All indices are base1
%
% NNO,TW,JD May 2010, updated Jan 2011
%
% See also SURFING_NODEIDXS2COORDS,SURFING_COORDS2LINVOXELIDXS,
% SURFING_CIRCLEROI, SURFING_SUBS2INDS, SURFING_INDS2SUBS

nverts=size(c1,2);

one_surf=numel(c2)==2;

dochecks=true; % do some input checking?
if dochecks
    if size(c1,1)~=3 || (~one_surf && size(c1,1)~=3) || size(f,1)~=3
        error('Coordinates and surfaces should be 3xP (or 3xQ)');
    end

    if ~isequal(unique(f(:))',1:nverts)
        error('Faces illegal; should contain all indices from 1 to %d',nverts)
    end
    
    if ~isequal(size(c1,2),size(c2,2)) && ~one_surf
        error('Coordinates C1 and C2 should have same number of elements');
    end
    
    mexfunctions={'surfing_uniqueidxsperrow','surfing_eucldist'};
    for k=1:numel(mexfunctions)
        whichmex=which(mexfunctions{k});
        [p,n,e]=fileparts(whichmex);
        if strcmpi(e,'.m')
            warning('SURFING:notoptimal',sprintf('%s was not compiled with mex.\nConsider running "mex %s.c" for increased speed.\n',n,n));
        end
    end
    
    if ~isstruct(voldef) || ~isfield(voldef,'mat') || ~isfield(voldef,'dim')
        error('voldef should be a struct with fields .mat and .dim');
    end
    
    if isfield(voldef,'mask') && prod(voldef.dim) ~= numel(voldef.mask)
        error('Volume mask provided, but number of elements does not match voldef.dim')
    end
end

if one_surf
    
    intermediatecoords=c1;
    
    start_pos=c2(1);
    stop_pos=c2(2);
    
    % make a copy
    c=c1;
    
    % compute normals
    c_normals=surfing_normals(c',f')';
    
    % compute offsets
    c1=bsxfun(@plus,c,start_pos*c_normals);
    c2=bsxfun(@plus,c,stop_pos*c_normals);
else
    
    % construct intermediate surface that is the average of the two surfaces
    intermediatecoords=squeeze(surfing_nodeidxs2coords(c1,c2,1:nverts,[1,0.5,0.5]));
end

    

usefixedradius=numel(circledef)==1;

if usefixedradius % number of voxels is not set
    circledef(2)=0;   
    radius=circledef(1);
elseif circledef(1)==0 || isnan(circledef(1)) % auto set radius in mm
    radius=10;
else 
    radius=circledef(1);
end

targetvoxcount=circledef(2); % target number of voxels (0 means no target, use radius only)


% set defaults 
if (nargin<6 || isempty(nodeidxs)),nodeidxs=1:nverts; end  % all nodes
if (nargin<7 || isempty(linedef)),linedef=[5 0 1];end     % 5 steps from surface 1 to surface 2
if (nargin<8 || isempty(distancemetric)),distancemetric='geodesic';end % geodesic distance metric
if (nargin<9 || isempty(progressstep)),progressstep=true; end       % show progress
if ~isnumeric(progressstep),progressstep=100;end % in steps of 100 nodes    

% number of nodes used as searchlight center
ncenters=numel(nodeidxs); 

% parameters in case targetvoxcount is set.
% strategy: use initial radius, then select all voxels with this radius.
% If too few voxels, increase radius, until enough voxels are selected
% now and then we update the value for the initial radius, so that most of the time 
% the initial radius is big enough but sometimes not. This is a tradeoff
% between having more nodes initially (with bigger radius; this slows down
% computing the geodesic distances), and increasing the radius sometimes
% which means re-computing geodesic distances
radiusgrow=1.5; % if radius is too small, multiply by this value and try again (see below)
updateradiuscount=floor(log2(ncenters));                                  % } Set how often we update the optimal radius, if targetvoxcount is set.
updateradiusat=repmat(2,1,updateradiuscount-3) .^ (4:updateradiuscount);  % } This is purely for faster execution. Update after 16, 32, ... nodes
updateradiusratio=0.80; % try to find the radius so that in 80% of the case we don't have to increase it. This is an emperical value.
updateradiusmax=1000; % just be be sure the radius does not grow infinitely

% NNO Sep 2011 added support for separate mask.centermask (the volume that 
% defines which voxels may contain center nodes) and mask.mask (the volume
% that defines which voxels can be selected). Also swapped the order of
% computing intermediate nodes and all nodes for obvious reasons.
% rationale: sometimes the mask may miss nodes that are at the edge of the
% volume

% find linear indices of voxels containing the center nodes. Values of NaN
% indicate outside the volume, and those are not used as center.
intermediatelinvoxidxs=surfing_coords2linvoxelidxs(intermediatecoords(:,nodeidxs),voldef); 
outsidemsk=isnan(intermediatelinvoxidxs);
if sum(outsidemsk)>0 && progressstep    
    warning('surfing_voxelselection:note','found %d / %d center nodes outside the volume, these will be ignored.', sum(outsidemsk), ncenters);
    %outsidemsk(:)=false;
    %NNO May 2011 still better to keep these nodes
    %NNO June 2011 or not!
    nodeidxs(outsidemsk)=NaN; %NNO Jan 2011
end
if isfield(voldef,'centermask')
    voldef=rmfield(voldef,'centermask');
end


% Find coordinates of points on or between the two surfaces.
% Lines between the surfaces are constructed according to LINEDEF.
% ALLCOORDS(I,J,K) is the I-th spatial coordinate (1,2,3)
% for the J-step along the line between the K-th node on surface 1 and 2
allcoords=surfing_nodeidxs2coords(c1,c2,1:nverts,linedef);

% Find the voxel indices corresponding to the coordinates of points above
% ALLLINVOXIDXS(I,K) contains the linear index for the voxel
% associated with node I for the K-th step. 
alllinvoxidxs=surfing_coords2linvoxelidxs(allcoords,voldef);

% For each row seperately, duplicates are replaced by zeros. 
unqlinvoxidxs=surfing_uniqueidxsperrow(alllinvoxidxs); 
clear alllinvoxidxs;

% set the node indices in random order, for better ETA estimation
nodeorder=randperm(ncenters);

% find the mapping from nodeidxs to the faces that contain the nodes
% this increases the speed of SURFING_SUBSURFACE dramatically
n2f=surfing_nodeidxs2faceidxs(f);

% construct mapping from linear to sub indices
lin2sub=surfing_inds2subs(voldef.dim,1:prod(voldef.dim));

% allocate space for the output
n2v=cell(1,ncenters);

% minimum and maximum values for voxel sub indices; or NaN if node is not
% selected
voxmin=NaN(ncenters,3);
voxmax=NaN(ncenters,3);

% number of voxels OR radius of the searchlights for each center node.
vORr=NaN(ncenters,1); % (number of) voxels OR radius chosen

ascenter=false(ncenters,1); %keep track which nodes were used as center

voxcountsum=0; % to keep track of average number of voxels

if isfield(voldef,'mask') && progressstep
    fprintf('Using %d / %d voxels in functional volume mask\n', sum(voldef.mask(:)~=0), numel(voldef.mask));
end
    
tic();
clock_start=clock();
prev_msg='';

for k=1:ncenters
    nodek=nodeorder(k);      %nodeorder ensures random order of visiting the nodes
    nodeidx=nodeidxs(nodek); %node index, in range 1:nverts

    % general while loop;
    % - if number of voxels is not given it exits after one iteration. 
    % - if number of voxels is given, then run one iteration, see if we
    %   found enough voxels, and if not we increase the radius and try
    %   again, until we selected enough voxels.

    radiusk=radius; % set initial radius
    ignorenode=outsidemsk(nodek);
    voxcount=0;
    while ~ignorenode % ignore node if outside the brain, or if radius has become too big
        
        % construct a circular ROI around node NODEIDX, and return the indices of
        % the nodes that are near to this node.
        if radius==0
            coordidxs=nodeidx;
            dist=0;
        else
            [coordidxs,dist]=surfing_circleROI(intermediatecoords,f,nodeidx,radiusk,distancemetric,n2f); 
        end

        % find voxel indices of voxels associated with the nodes in the circular ROI 
        
        linvoxidxs=unqlinvoxidxs(coordidxs,:);
        n2vk=unique(linvoxidxs(linvoxidxs>0));

        voxcount=numel(n2vk); % number of selected voxels
        if usefixedradius % select by radius; we're done, exit while-loop
            d=voxcount;
            break; % we're done for this node
        else
            
            % see how many voxels we have
            if voxcount<targetvoxcount 
                % not enough voxels; make radius bigger
                radiusk=radiusk*radiusgrow;
                
                if radiusk>updateradiusmax
                    ignorenode=true; % safety mechanism in case we cannot find
                    if progressstep
                        warning('surfing_voxelselection:note','could not find %d voxels for node %d, ignoring this node', targetvoxcount, nodeidx);
                    end
                end
                
                % we try again with this larger radius in the next iteration
            else 
                % we found enough voxels. 
                
                % sort the distance values to find the nodes that are
                % closest
                [foo,sidxs]=sort(dist);
                
                % linear indices associated with the nearest nodes, sorted
                % by distance
                linvoxidxs=unqlinvoxidxs(coordidxs(sidxs),:);
                
                % select approximately targetvoxcount voxels. nrows means
                % here the number of nodes associated with voxels that are
                % selected
                [n2vk,nrows]=surfing_selectkfirstidxs(targetvoxcount,linvoxidxs);
                
                % update voxelcount
                voxcount=numel(n2vk); 
                
                % distance of furthest node on the intermediate surface
                % (We don't take the (additional) distance
                % between a node on the center surface and the center of
                % the voxel into account.)
                d=dist(sidxs(nrows)); 
                break; % we're done for this node
            end
                    
        end
    end

    ignorenode=ignorenode || voxcount==0;
    
    if ~ignorenode
        % store results
        n2v{nodek}=uint32(n2vk(:))';
        ijk=lin2sub(n2vk,:);             % ijk indices (from linear indices)
        voxmin(nodek,:)=min(ijk,[],1); % minimum voxel indices
        voxmax(nodek,:)=max(ijk,[],1); % maximum voxel indices
        vORr(nodek)=d;    %maximum distance between node and center node
        ascenter(nodek)=true;    
    end
    
    
    % optionally show progress
    if progressstep
        voxcountsum=voxcountsum+voxcount;
        ascentercount=sum(ascenter);
        if k==1 || mod(k,abs(progressstep))==0 || k==ncenters; % show progress in beginning, every PROGRESSSTEP nodes, and at the end
            if progressstep<0, clc(); end
            tc=toc();
            eta=(ncenters-k)/k*tc;
            if usefixedradius
                rtxt=sprintf('r=%.2f, %.1f vox',radius,voxcountsum/ascentercount);
            else
                rtxt=sprintf('r=%.2f, %.1f vox',mean(vORr(ascenter)),voxcountsum/ascentercount);
            end
            msg=rtxt;
            prev_msg=surfing_timeremaining(clock_start,k/ncenters,msg,prev_msg);
            %fprintf('Completed %d / %d nodes (%.1f%%), %s, took %d sec, ETA %d sec\n',k,ncenters,k/ncenters*100,rtxt,round(tc),round(eta));
        end
    end
    
    % see if we should update the (hopefully close to optimal) radius, 
    % in case the target number of voxels is set.
    if ~usefixedradius && sum(updateradiusat==k) 
        ascentercount=sum(ascenter);
        radiusidx=round(updateradiusratio*ascentercount);
        if radiusidx>0
            dssort=sort(vORr(ascenter));
            radius=dssort(radiusidx); % find the radius in the updateradiusratio*100-th percentile 
            %if progressstep
            %    fprintf('After %d / %d nodes, radius set to %.2f (%d-th percentile)\n', k, ncenters, radius, round(updateradiusratio*100));
            %end
        end
    end
end

