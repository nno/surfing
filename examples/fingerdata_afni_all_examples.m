function fingerdata_afni_all_examples(just_one_surf,low_res_output)
% General illustration of 4 ways to do surface-based searchlight analysis
%
% fingerdata_afni_all_example(just_one_surf,low_res_output)
%
% This new example (as of May 2014) illustrates four ways to do searchlight
% analysis:
% - either using a low resolution surface to define the searchlight 
%   centers combined with a high resolution surface that delineate the 
%   grey matter; or using a high-resolution surface throughout. 
%   The advantage of a low-res surface is reduced execution time.
%   The disadantage is a loss in spatial resolution in the output map.
% - using a single surface (as provided by Caret or BrainVoyager) or
%   two surfaces (as provided by FreeSurfer)
%   Unlike the original fingerdata_afni_example_mvpa script, it uses the
%   output from pyref (see python/afni_anat_preproc.py), incorporates
%   the voxel selection step in this script, and uses merged hemispheres
%   (in which first nodes are stored from the left hemisphere, followed by
%   nodes in the right hemisphere). It also runs with just a single
%   hemisphere.
%
% The typical use case scenarios are:
% 1) BV/Caret users (just_one_surf=true)
%    (a) (low_res_output=false) use of single surface with same output 
%        resolution as the input resolution
%    (b) (low_res_output=true) use of single surface with lower output 
%        resolution as the input resolution. The output surface 
%        is downsampled using surfing_subsample_surface applied to the 
%        input surface 
% 2) FreeSurfer users (just_one_surf=false)
%    (a) use of pial and white surface with same output 
%        resolution as the input resolution
%    (b) (low_res_output=true) use of single surface with lower output 
%        resolution as the input resolution. The output surface 
%        is generated using MapIcosahedron (as is the input surface), with
%        the input surface having more linear divisions (and thus more 
%        nodes) than the output surface 
%
% Example:
%  >> fingerdata_afni_all_examples();
%
% NNO May 2014

if nargin==0
    % run all 4 examples in succession
    me=str2func(mfilename()); % make immune to renaming
    for just_one_surf=[true false]
        for low_res_output=[true false];
            me(just_one_surf,low_res_output);
        end
    end
    return
end

fprintf('\nStarting analysis; just_one_surf=%d, low_res_output=%d\n',...
                just_one_surf, low_res_output);
tic


% set parameters
radius=50; % searchlight size
radiusunit='vx'; % 'vx' (number of voxels) or 'mm' (circle radius)
hemi='m'; % merged hemispheres (left and right); alternatively, 'l' or 'r'

ld_high=64; % number of linear division for high-res surface

if low_res_output
    if just_one_surf
        downsample_niter=10; % Caret/BV case
    else
        ld_low=16; % FS case
    end
end

if just_one_surf
    % define offsets of points along normal vectors of nodes
    offset_inwards=-2; % start 2 mm 'inwards'
    offset_outwards=3; % stop 3 mm 'outwards'
    offsets=[offset_inwards offset_outwards];
end

% input directories
rootdir='/Users/nick/organized/_datasets/fingerdata-0.3/'; 
surfdir=[rootdir 'pyref/'];

% high res surfaces that define gray matter
% Note: if just_one_surface then a single surface is generated below
%       and these surface are not used anymore
pial_fn=sprintf('%s/ico%d_%sh.pial_al.asc', surfdir, ld_high, hemi);
white_fn=sprintf('%s/ico%d_%sh.smoothwm_al.asc', surfdir, ld_high, hemi);


in_dir=[rootdir 'glm/']; % input dir
out_dir=in_dir;           % output dir

beta_fn=[in_dir 'rall+orig'];

% load high-res surfaces
[v_pi, f_pi]=surfing_read(pial_fn);
[v_wh, f_wh]=surfing_read(white_fn);

% little sanity check
assert(isequal(f_wh,f_wh));

% define intermediate surface, if necessary
if just_one_surf || low_res_output
    v_im=(v_pi+v_wh)/2;
    f_im=f_wh;
end

    
if low_res_output    
    if just_one_surf
        out_infix=sprintf('%diter_%.1f-%.1f',downsample_niter,offsets);
        fprintf('Downsampling surface with %d iterations\n', ...
                                                downsample_niter);
        [v_src, f_src]=surfing_subsample_surface(v_im,f_im,...
                                            downsample_niter);
    else
        out_infix=sprintf('icolow%d',ld_low);
        
        % low res surface (from MapIcosahedron) that provides the 
        % the searchlight centers [FS case]
        source_fn=sprintf('%s/ico%d_%sh.intermediate_al.asc', ...
                    surfdir, ld_low, hemi);
        [v_src, f_src]=surfing_read(source_fn);
        fprintf('Source surface read from %s\n', source_fn);
    end
    
    % compute mapping from low to high res surface
    low2high=surfing_maplow2hires(v_src', v_im');
else
    if just_one_surf
        out_infix=sprintf('%.1f-%.1f',offsets);
    else
        out_infix='';
    end
    
    low2high=[]; % identity mapping
end
   
% sanity check: if just one surface remove pial and white from memory
% (this is the BV/Caret case)
if just_one_surf
    clear v_pi f_pi v_wh f_wh
end

% in the case of just one surface and a low-res surface, store new surfaces
% to disc
if just_one_surf && low_res_output
    source_fn=sprintf('%s/subs_ico%d_%dit_%sh.intermediate_al.asc',...
                out_dir, ld_high, downsample_niter, hemi);

    surfing_write(source_fn, v_src, f_src);
    fprintf('Output surface written to %s\n', source_fn);

    % also store an inflated version
    inflated_fn=sprintf('%s/ico%d_%sh.inflated_alCOMmedial.asc', ...
                                        surfdir, ld_high, hemi);

    
    [v_inf,unused]=surfing_read(inflated_fn);
    
    % apply mapping from low to high res to this surface
    v_src_inf=v_inf(low2high,:);
    
    % copy topology
    f_src_inf=f_src;
    
    
    source_inflated_fn=sprintf('%s/subs_ico%d_%dit_%sh.inflated_alCOMmedial.asc',...
                out_dir, ld_high, downsample_niter, hemi);
    surfing_write(source_inflated_fn, v_src_inf, f_src_inf);
    fprintf('Inflated surface written to %s\n', source_inflated_fn);

            
end
                                        

fn_out_data=sprintf('%s/%sh_ico%d_%s_%d%s',...
            out_dir,hemi,ld_high,out_infix,radius,radiusunit);

    

% define searchlight mapping
switch radiusunit
    case 'vx'
        circledef=[10 radius];
    case 'mm'
        circledef=radius;
    otherwise
        error('Unexpected unit for radius');
end

voldef=surfing_afni2spmvol(beta_fn);

% define voxel mask
opt=struct(); opt.Format='vector'; opt.Frames=1; % load first volume
[err,V,I]=BrikLoad(beta_fn,opt);
voldef.mask=V~=0; % only voxels in volume

% run voxel selection

fprintf('Starting voxel selection\n');
if just_one_surf
    n2v=surfing_voxelselection(v_im', offsets, f_im', circledef, voldef, low2high);
else
    n2v=surfing_voxelselection(v_pi', v_wh', f_pi', circledef, voldef, low2high);
end

fprintf('Loading functional data\n');
% prepare loading the EPI data
[b2s,unq,n2vs]=surfing_reducemapping(n2v); % find which voxels were ever selected
opt=struct(); 
opt.Voxels=unq; 
opt.Format='vector';
[err,V,I]=BrikLoad(beta_fn,opt); % load functional data

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

fprintf('Starting classification analysis\n');

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
S.node_indices=find(msk)-1; % base1->base0
S.data=r(msk,:);
S.stats={'Zscore()','None'};
S.labels={'Zacc','acc'};

surfing_write([fn_out_data '.niml.dset'],S);
toc
fprintf('Analyis completed\n\n');

