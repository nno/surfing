% This new example (as of Nov 2012) illustrates using a low resolution 
% surface to define the searchlight centers combined with and a high 
% resolution surface that delineate the grey matter.
% The advantage is that the script runs a lot faster
% The disadantage is a loss in spatial resolution.
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

ld_low=16;
ld_high=64;


% input directories
rootdir='/Users/nick/organized/_datasets/fingerdata-0.3/'; 
surfdir=[rootdir 'pyref/'];

% high res surfaces that define gray matter
pialfn=sprintf('%s/ico%d_%sh.pial_al.asc', surfdir, ld_high, hemi);
whitefn=sprintf('%s/ico%d_%sh.smoothwm_al.asc', surfdir, ld_high, hemi);

% low res surface that serves as searchlight center
sourcefn=sprintf('%s/ico%d_%sh.intermediate_al.asc', surfdir, ld_low, hemi);

indir=[rootdir 'glm/']; % input dir
outdir=indir;           % output dir

betafn=[indir 'rall+orig'];
fnout_data=sprintf('%s/%sh_ico%d-%d_cfy_%d%s',outdir,hemi,ld_low,ld_high,radius,radiusunit);

% load surfaces
[v_pi, f_pi]=surfing_read(pialfn);
[v_wh, f_wh]=surfing_read(whitefn);
[v_src, f_src]=surfing_read(sourcefn);

% define high-res intermediate surface
v_src_highres=.5 * (v_pi + v_wh); 

% find mapping from low to high resolution surface
low2high=surfing_maplow2hires(v_src', v_src_highres');

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
n2v=surfing_voxelselection(v_pi', v_wh', f_pi', circledef, voldef, low2high);

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

ncenters=numel(n2vs); % number of searchlights
rp=randperm(ncenters);
pred=zeros(1,nchunks*nclasses);
r=zeros(ncenters,2); % output result
msk=false(ncenters,1); % mask of nodes processed

clock_start=clock();
prev_msg='';
progress_step=100;

for k=1:ncenters
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
    
    if mod(k,progress_step)==0, 
        msg=sprintf('Mean z=%.3f, accuracy=%.3f ',mean(r(msk,:)));
        prev_msg=surfing_timeremaining(clock_start,k/ncenters,msg,prev_msg);
    end
end

S=struct();
S.node_indices=find(msk)-1;
S.data=r(msk,:);
S.stats={'Zscore()','None'};
S.labels={'Zacc','acc'};

surfing_write([fnout_data '.niml.dset'],S);





