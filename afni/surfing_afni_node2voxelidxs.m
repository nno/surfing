function idxs=surfing_afni_node2voxelidxs(Info,vertw,vertp,n,orientation)
% projects nodes on a surface to indices on the surface
%
% IDXS=SURFING_NODE2VOXELIDXS(INFO,WSURF,PSURF,N,ORIENT) takes an afni
% volume with info INFO and inner and outer surfaces WSURF and PSURF in
% orientation ORIENT, and returns N linear indices per node in IDXS. 
%
% WSURF and PSURF should be Px3 matrices for P nodes with x, y, z
% coordinates. The result IDX is a PxN matrix
%
% NNO Jan 2010

if nargin<4
    n=10;
end
if nargin<5
    orientation='LPI';
end

if nargin<3
    fprintf('Using test data!');
    [err,Info]=BrikInfo('/Users/psyuser/Documents/organized/105_observe_execute_actions2/sr906/anat+orig/sr906_frun01_preproc_raw+orig');
    nverts=20;
    m=20;
    vertw=repmat(rand(nverts,1)*m,1,3) + repmat([1 10 100],nverts,1);
    vertp=vertw;%1+repmat(rand(nverts,1),1,3) + repmat([1 10 100],nverts,1);
    n=8;
end

[nvertw,three1]=size(vertw);
[nvertp,three2]=size(vertp);

if nvertw ~= nvertp
    error('Number of vertices (%d and %d) do not match', nvertw, nvertp);
end
if three1 ~= 3 || three2 ~= 3
    error('Vertices should be specified in Px3 matrix');
end

% Relative position s is specified as going along the line from a wm (=0) to
% pial node (=1). 
steps_w2p=zeros(1,1,1,n);
steps_p2w=zeros(1,1,1,n);  

if n==1
    verts_weighted=0.5*(vertw+vertp);
%    steps_w2p(:)=0.5; %in between
%    steps_p2w(:)=0.5; %
else
    minpos=0;
    maxpos=1;
    relpos=minpos:(maxpos-minpos)/(n-1):maxpos;
    steps_w2p(:)=relpos;
    steps_p2w(:)=relpos(n:-1:1);


    verts=zeros(nvertw,3,2,'single'); %vertices from both surfaces
    verts(:,:,1)=vertw;
    verts(:,:,2)=vertp;
    verts_rep=repmat(verts,[1,1,1,n]);

    steps_rep=zeros(nvertw,3,2,n,'single');
    steps_rep(:,:,1,:)=repmat(steps_w2p,[nvertw,3,1,1]);
    steps_rep(:,:,2,:)=repmat(steps_p2w,[nvertw,3,1,1]);
    
    % weighted vertices, and reshaped in a Qx3 matrix
    verts_weighted=reshape(shiftdim(squeeze(sum(steps_rep .* verts_rep, 3)),2),nvertw*n,3);
end

% convert to voxel indices
[err,coordslin]=AFNI_XYZcontinuous2Index(verts_weighted,Info,orientation,1,true); %gets (i,j,k) voxel coordinates, 0-based
coordslin_base1=coordslin+1; % convert to 1-based indexing, as that's what matlab uses.
idxs=reshape(coordslin_base1,n,nvertw)';
idxs(isnan(idxs))=0;

