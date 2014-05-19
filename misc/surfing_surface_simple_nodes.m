function [pths,is_simple]=surfing_surface_simple_nodes(v,f)
% find 'simple nodes' in a surface
%
% [pths,is_simple]=surfing_surface_simple_nodes(v,f)
%
% Inputs:
%   v         Px3 coordinates for P vertices
%   f         Qx3 vertex indices for Q faces
%
% Outputs:
%   pths       PxM node indices, if each node is contained in at most M
%              faces. pths(K,:) contains as the first N elements a list of 
%              indices that form a path around node K, if such a path 
%              exists; the remaining entries are zero. 
%   is_simple  indicates whether a path exists around each node
%
% Notes:
%   - a path around a node K is a list of node indices I_1,...,I_N so that
%     the edges (I_1, I_2), (I_2, I_3), ... (I_{N-1}, I_N), (I_N, I_1)
%     are all contained in a face, that each of these faces contains 
%     node K, and that there are no other faces that contain node K.
%
% Example:
%  - consider the following surface, with arabic numerals indicating node
%    indices and roman numerals face indices:
%
%          2 ----- 7
%         / \     / \
%        / I \ VI/V  \
%       /     \ /     \
%      3 ----- 1 ----- 6
%       \     / \     /
%        \ II/III\ IV/
%         \ /     \ /
%          4 ------5
%
%    this surface can be described by
%    >> v=[0 -1 -2 -1 1 2 1;0 -2 0 2 2 0 -2;0 0 0 0 0 0 0]';
%    >> f=[1 1 1 1 1 1; 2 3 4 5 6 7; 3 4 5 6 7 2]';
%
%    In this surface, only node 1 is a simple node as there is a path
%    around it 2-3-4-5-6-7(-2) that covers all faces (I-VI) that contain 
%    node 1 while the path starts and ends with the same node.
%    All other nodes are not simple because they don't have such a path.
%
%    >> [pths,is_simple]=surfing_surface_simple_nodes(v,f)
%    pths =
% 
%      2     3     4     5     6     7
%      0     0     0     0     0     0
%      0     0     0     0     0     0
%      0     0     0     0     0     0
%      0     0     0     0     0     0
%      0     0     0     0     0     0
%      0     0     0     0     0     0
% 
% 
%    is_simple =
% 
%      1
%      0
%      0
%      0
%      0
%      0
%      0

nv=size(v,1);
nf=size(f,1);

% mapping from nodes to faces
v2f=surfing_nodeidxs2faceidxs(f');

% maximum number of faces per node
nmax=size(v2f,2);

% space for output
pths=zeros(nv,nmax);
is_simple=false(nv,1);

% allocate temporary space for path
pth=zeros(1,nmax);

% keep track of path lengths
npths=zeros(nv,1);

% keep track how often each face was visited
col_counter=zeros(1,nmax);

clock_start=clock();
prev_msg='';
progress_step=5000;


for k=1:nv
    % for each node get the face indices that contain it
    fs=v2f(k,:);

    % node indices (Px3) array of the faces that contain node k
    ft=f(fs(fs>0),:)';

    n=size(ft,2);
    ns=(1:n)';

    a=find(ft==k);
    fa=ft(mod(a,3)+1+ns*3-3)';  % successor   } of each node, assuming
    fb=ft(mod(a+1,3)+1+ns*3-3); % predecessor } there is a circular path

    col=1;
    col_counter(:)=0; % reset counter
    
    has_path=true;
    
    for j=1:n
        pth(j)=fa(col);

        col=find(fb(col)==fa);

        if numel(col)~=1
            has_path=false;
            break;
        end
        
        col_counter(col)=col_counter(col)+1;
    end

    if has_path && all(col_counter(ns)==1)
        pths(k,ns)=pth(ns);
        is_simple(k)=true;
        npths(k)=n;
    end  
    
    if mod(k,progress_step)==0 || k==nv
        s=is_simple(1:k);
        msg=sprintf('%.1f%% simple nodes', sum(s)*100/nv);
        prev_msg=surfing_timeremaining(clock_start,k/nv,msg,prev_msg); 
    end
end

% special case: neighbors of nodes with three nodes around them
% these nodes are not simple
% example: surfaces with faces 1-4-2, 1-2-5, 2-4-3, 2-3-5, 3-4-5
% 

triples_msk=npths==3;
pth2node=surfing_invertmapping(pths);

non_simple=sum(pths>0,2)<3;

for col=1:3
    % get pivot nodes
    pivots=pths(triples_msk,col);

    % get paths of pivot nodes
    p=pths(pivots,:);
    
    % get nodes on the path before and after the pivot
    a=bsxfun(@eq,pths(triples_msk,mod(col,3)+1),p);
    b=bsxfun(@eq,pths(triples_msk,mod(col+1,3)+1),p);
    
    % predecessor and sucessor should not be equal to pivot
    m=any(a,2) & any(b,2);
    
    non_simple(pivots(m))=true;
end

pths(non_simple,:)=0;
is_simple(non_simple)=false;
        