function [vertices,faces]=surfing_generate_planar_surface(nx, ny, origin, x1, y1)
% generate a surface in the shape of a plane
%
% [vertices,faces]=surfing_generate_planar_surface(nx, ny, origin, x1, y1)
%
% Inputs:
%   nx          number of vertices in first dimension
%   ny          number of vertices in second dimension
%   origin      (optional) 1x3 vector with coordinates of origin. 
%               Default: [0,0,0]
%   x1          (optional) 1x3 vector with coordinates of second vertex
%               (the origin being the first) in the first dimension. 
%               Default: [1,0,0]
%   y1          (optional) 1x3 vector with coordinates of second vertex
%               (the origin being the first) in the second dimension
%               Default: [0,1,0]
%
% Output:
%   vertices    Px3 coordinates, with P = nx * ny
%   faces       Qx3 faces, with Q=2*(nx-1)*(ny-1)
%
% Example:
%     % generate a surface 
%     [v,f]=surfing_generate_planar_surface(5,2,[0,0,0],[1,0,0],[0,1,0])
%     > v =
%     > 
%     >      0     0     0
%     >      0     1     0
%     >      1     0     0
%     >      1     1     0
%     >      2     0     0
%     >      2     1     0
%     >      3     0     0
%     >      3     1     0
%     >      4     0     0
%     >      4     1     0
%     > 
%     > 
%     > f =
%     > 
%     >      1     2     3
%     >      4     3     2
%     >      3     4     5
%     >      6     5     4
%     >      5     6     7
%     >      8     7     6
%     >      7     8     9
%     >     10     9     8
%
% Notes:
%   - this function generates a surface consisting of rectangles, with
%     three corners of the first rectangle being the origin, x1 and y1.
%     Each rectangle consists of two triangles
%
% NNO Sep 2015

if nargin<3
    origin=[0,0,0];
end

if nargin<4
    x1=[1,0,0];
end

if nargin<5
    y1=[0,1,0];
end


cellfun(@ensure_is_three_vec,{origin,x1,y1});

vertices=zeros(nx*ny,3);
faces=zeros(2*(nx-1)*(ny-1),3);

for i=1:nx
    for j=1:ny
        vpos=(i-1)*ny+j;
        vertices(vpos,:)=origin+(i-1)*x1+(j-1)*y1;
        if i<nx && j<ny
            p=vpos;
            q=vpos+1;
            r=vpos+ny;
            s=vpos+ny+1;
            
            fpos=((i-1)*(ny-1)+j)*2-1;
            faces(fpos,:)=[p,q,r];
            faces(fpos+1,:)=[s,r,q];
        end
    end
end

function ensure_is_three_vec(x)
    if ~(isnumeric(x) && ...
            numel(x)==3 && ...
            numel(size(x))==2)
        error('vectors must have three numeric elements');
    end
    