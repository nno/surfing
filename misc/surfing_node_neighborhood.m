function [node_pth,face_pth,theta]=surfing_node_neighborhood(v,f,pivot,v2f)
% determine neighborhood information for a node
%
% [node_pth,face_pth,is_convex]=surfing_node_neighborhood(v,f,pivot,v2f)
%
% Inputs:
%  v          Px3 coordinates for P nodes
%  f          Qx3 node indices for Q faces
%  pivot      node index for which neighborhood information is returned
%  v2f        PxM mapping from nodes to faces, if each node is contained in
%             at most M faces (optional). It should be the output from
%             surfing_invertmapping(f) and can be provided as argument to
%             avoid recomputing the mapping upon every call of this
%             function. If not given it is assigned the output from
%             surfing_invertmapping(f)
%
% Outputs:
%  node_pth   Kx1 node indices that form a path around the neighbors of the
%             pivot node. Thus, if node_pth=[N_1,...N_K) and N_0==N_K, then
%             N_* is the set of neighbors of the pivot node and all pairs
%             [N_(J-1), N_J)] are an edge on the surface. If no such path
%             exists it is set to the empty list [].
%  face_pth   Kx1 face indices that contain the pivot node
%  is_convex  boolean indicating whether the surface surrounded by the
%             nodes in node_pth is convex. If no path exists in node_pth
%             this value is always false
%
% NNO Apr 2014


% compute mapping from vertices to faces if not given
if nargin<4 || isempty(v2f)
    v2f=surfing_invertmapping(f);
end


% triangles around pivot node
tris=v2f(pivot,:);
tris=tris(tris>0);
ntri=numel(tris);  % number of triangles

% nodes in triangles
fp=f(tris,:);

% allocate space
node_pth=zeros(ntri,1);  % path of nodes around pivot
face_pth=zeros(ntri,1);  % path of faces around pivot
pth_pos=0;

visited=false(ntri,1);
has_no_path=false;


% row and column in fp node array
row=find(~visited,1);
col=mod(find(fp(row,:)==pivot),3)+1; % one next after pivot

% find path around pivot that visits all its neighbors in order
while ~visited(row)
    visited(row)=true;
    rowval=fp(row,:);
    pval=rowval(col);

    % store path information
    pth_pos=pth_pos+1;
    node_pth(pth_pos)=pval;
    face_pth(pth_pos)=tris(row);

    % take remaining element
    % in case the surface is poorly behaved (i.e. inconsistent orientation)
    % find the 'other' element rather than the next one
    %qval=fp(row,mod(col,3)+1); % next node element
    qval=rowval(rowval~=pivot & rowval~=pval);

    msk=fp==qval;
    msk(row,:)=false;

    % if not a closed surface, there may be no path
    if sum(msk(:))~=1
        visited(:)=false;
        break;
    end
    % determine position of next node on path
    col=find(sum(msk,1));
    row=find(sum(msk,2));
end

if ~all(visited)
    % maybe use this in the future
    has_no_path=true;
end


if has_no_path
    node_path=[];
    face_path=[];
    theta=NaN(ntri,1);
else
    is_convex=true;
    theta=zeros(ntri,1);

    vec_as=zeros(ntri*2,3);
    vec_bs=zeros(ntri*2,3);

    for k=1:ntri
        idxs=node_pth(mod(k-1+(0:2),ntri)+1);
        % points: pivot p; neighbors point a, b and c on the path
        % ensure angle(a,b,p)+angle(c,b,p)<pi
        xyz=v([pivot; idxs],:);

        bp=xyz(1,:)-xyz(3,:);
        ab=xyz(2,:)-xyz(3,:);
        bc=xyz(4,:)-xyz(3,:);

        % two internal angles
        %t1=angle(ab,bp);
        %t2=angle(bc,bp);
        %theta(k)=abs(t1)+abs(t2);

        %ts=angle([ab; bc],[bp; bp]);
        %theta(k)=sum(abs(ts));
        idxs=(k-1)*2+(1:2);
        vec_as(idxs,:)=[ab; bc];
        vec_bs(idxs,:)=[bp; bp];
    end

    thetas=abs(angle(vec_as, vec_bs));
    theta=sum(reshape(thetas,2,[]),1)';
end

function theta=angle(a,b)
    % helper function to compute angle between two vectors

    %theta=acos((a*b')/sqrt((a*a')*(b*b')));
    theta=acos( sum(a.*b,2)./sqrt (sum(a.^2,2).*sum(b.^2,2)));
