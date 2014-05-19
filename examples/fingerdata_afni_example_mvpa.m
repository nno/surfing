radius=50;
radiusunit='vx'; % vx or mm
hemi='l'; % left hemisphere                    
centercount=NaN; % number of searchlight centers; NaN for all

% input directories
rootdir='/Users/nick/organized/_datasets/fingerdata-0.3/'; %'../../surfing_exampledata/';
indir=[rootdir 'glm/']; % input dir
outdir=indir;

voxselfn=sprintf('%s/%sh_voxsel_%d%s',indir,hemi,radius,radiusunit);
betafn=[indir 'rall+orig'];
fnout=sprintf('%s/%sh_cfy_%d%s',outdir,hemi,radius,radiusunit);

% load voxel selection data
R=load(voxselfn);

[b2s,unq,n2vs]=surfing_reducemapping(R.n2v); % find which voxels were ever selected
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

S=struct();
S.node_indices=R.centernodeidxs(msk)-1; % base0
S.data=r(msk,:);
S.stats={'Zscore()','None'};
S.labels={'Zacc','acc'};

afni_niml_writesimple([fnout '.niml.dset'],S);





