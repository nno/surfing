% Example function to demonstrate surface-based voxel selection.
% 
% This function uses two surfaces (pial and white matter, from Freesurfer).
% On the surface, a "center" node is selected randomly, nodes nearby this
% center node are selected, the corresponding voxels are found, and a value
% of one is added to these voxels. This is repeated many (default: 1000)
% times, resulting in a volume where (almost or all) voxels falling in the 
% left cortex have a value of one or more, but nodes outside the cortex 
% a value of zero. 
%
% This function requires some example data: a volume and two surfaces
%
% NNO May 2010

ncenters=1000; % how many searchlights
circledef=[10 200]; % initial searchlight radius r=10, find 200 voxels per searchlight.
                    % use [10 0] for a radius of 10, with variable number of voxels.

% input directories
rootdir='/home/tobias/Desktop/latest_surfing/surfing_exampledata/'; %'../../surfing_exampledata/';
fsdir=rootdir; % freesurfer dir
bdir=rootdir; %  beta dir cortex

% input files
fns1=[fsdir 'ico100_lh.pial_al.asc'];     % freesurfer pial surface (ascii format) 
fns2=[fsdir 'ico100_lh.smoothwm_al.asc']; % freesurfer white surface (ascii format) 
% (to read binary freesurfer surfaces, use freesurfer_read_surf)

betatype='analyze';
switch betatype
    case 'analyze'
        fnb=[bdir 'epi.nii'];               
        V=spm_vol(fnb);
    case 'afni'
        fnb=[bdir 'epi+orig'];               
        [err,I]=BrikInfo(fnb);
        V=surfing_afni2spmvol(I);
    otherwise
        error('unknown betatype %s', betatype);
end


% output 
outdir=rootdir;
fnout=[outdir 'selected'];

% (Caret based, not used)
% %----go to the glm first level directory and save the images there
% cd(N.glmDIR);
% %----load the surface mask
% s1=caret_load(N.Right.coordImage_pial);
% s2=caret_load(N.Right.coordImage_white);
% f1=caret_load(N.Right.topoFile);                          

% load freesurfer ASCII surfaces
[c1,f]=freesurfer_asc_load(fns1);
[c2,f_]=freesurfer_asc_load(fns2);

% transpose for voxelselection
c1=c1';
c2=c2';
f=f';

nverts=size(c1,2); % number of vertices
allnodeidxs=randperm(nverts); % random permutation
nodeidxs=allnodeidxs(1:ncenters); % random selection of center indices
voldef=V;

% optional arguents (these values are the defaults, but shown for
% instructory purpose)
linedef=[5 0 1]; % lines from surface 1 to surface 2: 5 steps, start at surface 1, end at surface 2 
distancemetric='geodesic';               
progressstep=100; % show progress every 100 voxels; use 0 for no progress report

fprintf('Starting voxel selection, %d center nodes\n',ncenters);
[n2v,mn,mx,ds]=surfing_voxelselection(c1,c2,f,circledef,voldef,nodeidxs,linedef,distancemetric,progressstep);
fprintf('Completed voxel selection, %d center nodes\n',ncenters);

% output:
%  n2v     contains the node to linear voxel indices for NODEIDXS
%  mn, mx  the minimum and maximum sub voxel indices for each searchlight
%  ds      searchlight radius (contant nr of voxels), or number of voxels (constant radius) 


% prepare for storing volume data
K=zeros(V.dim); % number of times each voxel was selected by different searchlights

% run over all node indices
% more efficient (but less instructive) would be to use nodeidxs
for k=1:numel(n2v)
    linidxs=n2v{k}; % contains a 3xP_k matrix, if for node k there were P_k nodes selected
    if numel(linidxs)==0
        continue;
    end
    % use linear indices
    
    K(linidxs)=K(linidxs)+1;
end

% write results to volume
switch betatype
    case 'analyze';
        V.fname=[fnout '.nii'];
        spm_write_vol(V,K);
    case 'afni'
        opt=struct();
        opt.Prefix=[fnout '+orig'];
        opt.OverWrite='y';
        I.DATASET_RANK(2)=1;
        I.BRICK_LABS='selected';
        I=rmfield(I,'BRICK_STATS');
        I=rmfield(I,'BRICK_TYPES');
        I=rmfield(I,'BRICK_FLOAT_FACS');
        [e,E,I2]=WriteBrik(K,I,opt);
end
        
        