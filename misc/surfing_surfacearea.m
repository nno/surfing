function [node2area,face2area]=surfing_surfacearea(c,f,n2f)
% Computes node and face areas of a surface
%
% [NAREA,FAREA]=SURFING_SURFACEAREA(C,F[N2F])
% INPUTS:
%   C         Px3 (or 3xP) coordinates for P nodes
%   F         Qx3 (or 3xQ) node indices for Q faces
%   N2F       Optional node to face mapping (from SURFING_NODEIDXS2FACEIDXS)
%             If omitted, this mapping is computed on the fly.
% OUTPUTS:
%   NAREA     Px1 node area (based on faces containing the node)
%   FAREA     Qx3 face area
%
% NNO Oct 2010

if size(c,2) ~= 3, c=c'; end
if size(f,2) ~= 3, f=f'; end

[nv,three]=size(c);
[nf,three_]=size(c);

if three_ ~= 3 || three ~= 3, error('Coordinates and faces should be 3xP and 3xQ'); end

% coordinates of triangle nodes
p=c(f(:,1),:);
q=c(f(:,2),:);
r=c(f(:,3),:);

% vectors of two sides of triangles
pq=p-q;
pr=p-r;

% face area: .5*sqrt(   |a|^2*|b|^2 - (a.b)^2    ) 
% where a=pq, b=pr and (a.b) means dot product
face2area=.5*sqrt(sum(pq.^2,2).*sum(pr.^2,2)-sum(pq.*pr,2).^2);

if nargin<3
    % find node to face maping
    n2f=surfing_invertmapping(f);
end

maxn2f=size(n2f,2);
node2area=zeros(nv,1);
for k=1:maxn2f
    n2fk=n2f(:,k); % not all nodes may be in the same number of faces
    msk=n2fk>0 & ~isnan(n2fk); 
    node2area(msk)=node2area(msk)+face2area(n2fk(msk))/3; % even split between three nodes in the triangle
end
    
