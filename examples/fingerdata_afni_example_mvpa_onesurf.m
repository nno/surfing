% This example (as of Nov 2012) illustrates using a a single
% high-resolution surface. 
% A use case is surfaces from caret or brainvoyager.
% In this example the single surface is the intermediate surface 
% from freesurfer (i.e. based on the node-wise average of the pial and 
% white surfaces), extended 2mm 'inwards' (towards the white matter) and 
% 3mm 'outwards' (towards the skull). Note that these values are quite
% arbitrary and are not necessarily optimal (based on whatever criterion).
%
% The extensions are made using the node-wise normal vectors; it may be the
% case the for brainvoyager the values must be flipped as it uses a
% different handednesss than freesurer.
%
% Unlike the original fingerdata_afni_example_mvpa script, it uses the
% output from pyref (see python/afni_anat_preproc.py), incorporates
% the voxel selection step in this script, and uses merged hemispheres
% (in which first nodes are stored from the left hemisphere, followed by
% nodes in the right hemisphere).
%
% NNO Nov 2012

radius=50;
radiusunit='vx'; % vx or mm
hemi='m'; % merged hemispheres (left and right)
centercount=NaN; % number of searchlight centers; NaN for all

% use a intermediate resolution surface
% for publication-quality results, 64 or 128 are more suitable
ld_high=64;

% define offsets of points along normal vectors of nodes
% BV users: inwards and outwards may have to be swapped due to
% handedness-differences between FS and BV.
offset_inwards=-2; % start 2 mm 'inwards'
offset_outwards=3; % stop 3 mm 'outwards'
offsets=[offset_inwards offset_outwards];

% input directories
rootdir='/Users/nick/organized/_datasets/fingerdata-0.3/'; 
surfdir=[rootdir 'pyref/'];

% high res surfaces that define gray matter
pialfn=sprintf('%s/ico%d_%sh.pial_al.asc', surfdir, ld_high, hemi);
whitefn=sprintf('%s/ico%d_%sh.smoothwm_al.asc', surfdir, ld_high, hemi);

indir=[rootdir 'glm/']; % input dir
outdir=[rootdir 'output/'];

betafn=[indir 'rall+orig'];
fnout=sprintf('%s/%sh_ico%d_single-surface_cfy_%d%s',...
                        outdir,hemi,ld_high,radius,radiusunit);

% load surfaces
[v_pi, f_pi]=freesurfer_asc_load(pialfn);
[v_wh, f_wh]=freesurfer_asc_load(whitefn);

% define high-res intermediate surface
% Note: with Caret and BV there is just one surface, which would be 'v' as
% defined here; in this example it is pretended there are no white or pial
% surfaces.
v=.5 * (v_pi + v_wh); 

% take the topology from the pial surface (should be identical to white
% surface)
f=f_pi; 

% in this example: ensure we're not using the pial and white surfaces 
% 'by accident'
clear v_pi f_pi v_wh f_wh

% define searchlight mapping
if ~strcmp(radiusunit,'vx')
    error('Unexpected unit for radius');
end
circledef=[10 radius];
voldef=surfing_afni2spmvol(betafn);

% define voxel mask
opt=struct(); opt.Format='vector'; opt.Frames=1; % load first volume
[err,V,I]=BrikLoad(betafn,opt);
voldef.mask=V~=0; % only voxels in volume

% run voxel selection
% note: the second argument is the offsets
n2v=surfing_voxelselection(v', offsets, f', circledef, voldef);

% prepare loading the EPI data
[b2s,unq,n2vs]=surfing_reducemapping(n2v); % find which voxels were ever selected
opt=struct(); 
opt.Voxels=unq; 
opt.Format='vector';
[err,V,I]=BrikLoad(betafn,opt); % load functional data

% set up simple 2 class cross validation
nchunks=16;
nclasses=2;

% classes (index and middle finger press), alternating subbriks
classes=repmat(1:nclasses,1,nchunks); 

tridxs=cell(nchunks,1); 
teidxs=cell(nchunks,1);

for j=1:nchunks
    teidxs{j}=(j-1)*nclasses+(1:nclasses);
    tridxs{j}=setdiff(1:(nclasses*nchunks),teidxs{j});
end

% standard deviation of binomial distribution
binosd=sqrt((nclasses-1)/nclasses^2/(nchunks*nclasses));

n=numel(n2vs); % number of searchlights
rp=randperm(n);
pred=zeros(1,nchunks*nclasses);
r=zeros(n,2); % output result
msk=false(n,1); % mask of nodes processed
tic();
for k=1:n
    centeridx=rp(k);
    n2vk=n2vs{centeridx};
    
    if numel(n2vk)<10
        continue;
    end
    
    Vk=V(n2vk,:); % data in k-th searchlight
    
    pred(:)=0;
    for j=1:nchunks
        pred(teidxs{j})=classify_lda_KclassesQuicker(...
            Vk(:,tridxs{j}),classes(tridxs{j}),Vk(:,teidxs{j}));
    end
    
    acc=sum(pred==classes)/numel(classes); % accuracy
    zacc=(acc-1/nclasses)/binosd; % z score
        
    r(centeridx,:)=[zacc acc];
    msk(centeridx)=true;
    
    if mod(k,100)==0, 
        msg=sprintf('Mean z=%.3f, accuracy=%.3f ',mean(r(msk,:)));
        surfing_timeremaining(k/n,msg); 
    end
end

% store surface results
S=struct();
S.node_indices=find(msk)-1;
S.data=r(msk,:);
S.stats={'Zscore()','None'};
S.labels={'Zacc','acc'};

afni_niml_writesimple([fnout '.niml.dset'],S);

% Volume output: how often each voxel was in a searchlight
K=zeros(voldef.dim); 

% run over all node indices
for k=1:n
    linidxs=n2v{k};  
    K(linidxs)=K(linidxs)+1;
end

% write brain mask
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



