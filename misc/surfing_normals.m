function [node_normals, face_normals]=surfing_normals(v,f)
% compute normals of nodes and faces
%
% Inputs:
%   v              Px3 node coordinates
%   f              Qx3 face indices
%
% Returns:
%   node_normals   Px3 node normals
%   face_normals   Qx3 face normals
%
% NNO Feb 2014


% compute face normals
a=v(f(:,1),:);
b=v(f(:,2),:);
c=v(f(:,3),:);

ab=a-b;
ac=a-c;

abXac=cross(ab,ac);

face_normals=bsxfun(@rdivide,abXac,surfing_eucldist(abXac',[0 0 0]'));

% compute node normals
nv=size(v,1);
n_sum=zeros(nv,3);

i=surfing_invertmapping(f);
n_max=size(i,1);
for k=1:3
    msk=i(:,k)>0;
    n_sum(msk,:)=n_sum(msk,:)+face_normals(i(msk,k),:);
end
    
node_normals=bsxfun(@rdivide,n_sum,surfing_eucldist(n_sum',[0 0 0]'));
