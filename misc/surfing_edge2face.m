function edge2face=surfing_edge2face(v,f)
% determine mapping from edges to faces
%
% edge2face=surfing_edge2face(v,f)
%
% Input
%   f              Qx3 face indices
%
% Output:
%   edge2face      PxP sparse matrix (with P=max(f(:)) the number of nodes)
%                  so that edge2face(i,j)==k means that the edge from i to
%                  j contains face k. Note that edge2face(j,i) contains the
%                  other face shared with the edge.
%
% Notes:
%   - if the surface has no consistent orientation, or an edge is present
%     in multiple faces, an error is thrown
%
% NNO Apr 2014

nf=size(f,1);
nv=max(f(:));

empty=@()zeros(nf*3,1);
ii=empty(); % row
jj=empty(); % column
face=empty(); % face indices
face_count=zeros(nf*3,1); % number of faces for each edge

for col=1:3
    % fill in edges from columns col and (col+1) mod 3
    row=(col-1)*nf+(1:nf);
    ii(row)=f(:,col);
    jj(row)=f(:,mod(col+1,3)+1);
    face(row)=1:nf;
    face_count(row)=face_count(row)+1;
end

% number of faces associated with each edge
% if more than one face for an edge an error is thrown
edge2face_count=sparse(ii,jj,face_count,nv,nv,nf*3);
mx=max(max(edge2face_count(:)))+0;

if mx>1
    [i,j]=find(edge2face_count>1,1);

    % find offending rows
    msk=ii==i & jj==j;
    assert(sum(msk)>=2);

    error('duplicate edge (%d,%d) in faces %d and %d',...
                    i,j,face(find(msk,2)));

end

edge2face=sparse(ii,jj,face,nv,nv,nf*3);
