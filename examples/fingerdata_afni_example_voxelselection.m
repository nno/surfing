% Example function to demonstrate surface-based voxel selection.
% 
% This function uses two surfaces (pial and white matter, from Freesurfer)
% to define circular regions of approximately 50 voxels using a geodesic
% distance measure.
%
% This function is aimed at the ''fingerdata' dataset
%
% NNO Feb 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radius=50;
radiusunit='vx'; % vx or mm
hemi='l'; % left hemisphere                    
centercount=NaN; % number of searchlight centers; NaN for all

% input directories
%rootdir='/Volumes/organized/116_surfing_toolbox/fingerdata-0.2/'; 
rootdir='/Users/nick/organized/_datasets/fingerdata-0.3/';
fsdir=[rootdir 'ref/']; % freesurfer dir
bdir=[rootdir 'glm/']; %  beta dir cortex
outdir=bdir; % output dir

% input files (ASCII freesurfer)
fns1=[fsdir 'ico100_' hemi 'h.pial_al.asc'];     % freesurfer pial surface (ascii format) 
fns2=[fsdir 'ico100_' hemi 'h.smoothwm_al.asc']; % freesurfer white surface (ascii format) 
% (to read binary freesurfer surfaces, use freesurfer_read_surf)

fnb=[bdir 'epiref+orig'];    
[err,I]=BrikInfo(fnb);
voldef=surfing_afni2spmvol(I); % make SPM-like struct

opt=struct(); opt.Format='vector'; opt.Frames=1; % load first volume
[err,V,I]=BrikLoad(fnb);
voldef.mask=V>0; % only voxels in volume

switch radiusunit
    case 'vx'
        circledef=[10 radius]; % initial radius 10mm
    case 'mm'
        circledef=radius;
    otherwise
        error('Illegal radius unit %s', radiusunit);
end

% output 
fnout=sprintf('%s/%sh_voxsel_%d%s',outdir,hemi,radius,radiusunit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load freesurfer ASCII surfaces
[c1,f]=freesurfer_asc_load(fns1);
[c2,f_]=freesurfer_asc_load(fns2);

% transpose for voxelselection
c1=c1';
c2=c2';
f=f';

nverts=size(c1,2); % number of vertices in a surface

if isempty(centercount) || isnan(centercount)
    centernodeidxs=1:nverts; % use all nodes as center
    centercount=nverts;
else
    allnodeidxs=randperm(nverts); % random permutation
    % random selection of center indices
    centernodeidxs=allnodeidxs(1:min(nverts,centercount)); 
    centercount=numel(centernodeidxs);
    fnout=sprintf('%s%dnodes',fnout,centercount);
end

% optional arguents 
% (these values are the defaults, but shown for instructory purpose)

% lines from surface 1 to surface 2: 5 steps, start at surface 1, end at
% surface 2 
linedef=[5 0 1]; 
distancemetric='geodesic'; 
progressstep=100; % show progress every 100 voxels; use 0 for no progress report


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do voxel selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Starting voxel selection, %d center nodes\n',centercount);
[n2v,mn,mx,ds]=surfing_voxelselection(c1,c2,f,circledef,voldef,...
               centernodeidxs,linedef,distancemetric,progressstep);
fprintf('Completed voxel selection, %d center nodes\n',centercount);

% output:
%  n2v     contains the node to linear voxel indices for NODEIDXS
%  mn, mx  the minimum and maximum sub voxel indices for each searchlight
%  ds      searchlight radius or number of voxels 

% save voxel selection output
save(fnout,'n2v','mn','mx','ds','centernodeidxs','radius','radiusunit');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More output for instructive purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save how many voxels were selected for each node for both (1) each 
% center node on the surface, and (2) each selected voxel in the volume

% Surface output: how many voxels in each searchlight
scount=zeros(centercount,1);
for k=1:centercount
    scount(k)=numel(n2v{k});
end

S=struct();
S.data=[scount+.0001,ds]; % ensure these are saved as floats 
S.data(isnan(S.data))=0; % SUMA does not appreciate NaNs
S.node_indices=centernodeidxs-1;
S.stats={'None','None'};
S.labels={'Nvoxels','radius'};
afni_niml_writesimple(S,[fnout '.niml.dset']); % write results

% (2) Volume output: how often each voxel was in a searchlight
K=zeros(voldef.dim); 

% run over all node indices
for k=1:centercount
    linidxs=n2v{k};  
    K(linidxs)=K(linidxs)+1;
end

% write results 
opt=struct();
opt.Prefix=[fnout '+orig'];
opt.OverWrite='y';

% take volume definition I from input brik
I.DATASET_RANK(2)=1;
I.BRICK_LABS='selected';
I=rmfield(I,'BRICK_STATS');
I=rmfield(I,'BRICK_TYPES');
I=rmfield(I,'BRICK_FLOAT_FACS');
[e,E,I2]=WriteBrik(K,I,opt);

        
        