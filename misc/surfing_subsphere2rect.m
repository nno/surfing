function [idxs,nx,ny]=surfing_subsphere2rect(v,io,ix,iy,nx,ny,progressstep)
% Projects part of a spherical surface to a rectangle
%
% [IDXS,NX,NY]=SURFING_SUBSPHERE2RECT(V,IO,IX,IY,NX,NY)
% INPUTS:
%   V           Px3 surface coordinates (x,y,z) for P nodes
%   IO,IX,IY    Node indices (base 1) for origin, x and y edges of
%               rectangle. If ABCD is a rectangle with sides AB, BC, CD,
%               and AC, then A=V(IO,:), B=V(IX,:), and D=V(IY,:).
%   NX,NY       Number of points on the sides of the rectangle 
%               If NY is omitted, is is computed to give a similar step
%               size in the Y direction as the X direction
% OUTPUTS:
%   IDXS        NXxNY matrix with linear indices of the neareast node of
%               the surface relative to the rectangle
%   NX,NY       Size of rectangle
%
% Method: Construct a rectangle based on IO,IX,IY, subdivide the sides,
% and compute the nearest points between those on the rectangle and the
% surface. Subdivision is currently based on distance, not angle. 
% Currently it is assumed that the input surface is an AFNI ascii spherical
% surface; the program may crash if those assumptions are not met.
%
% NNO Dec 2010

if nargin<6,ny=[];end
if nargin<7, progressstep=1e3; end

% If V is a filename of an freesurfer ASCII surface, try to load it
if ischar(v) && exist(v,'file') && strcmp(v(end+(-3:0)),'.asc')
    [v,f]=freesurfer_asc_load(v);
end

% check size of input
[s1,s2]=size(v);
if s1==3
    v=v';
elseif s2~=3
    error('V should be Px3')
end

v=bsxfun(@minus,v,mean(v));

% check that surface is a sphere with radius 100
dd=sum(v.^2,2);
if abs(dd-1e4)>1
    error('Expected sphere with radius 100');
elseif max(mean(v))>1e-3
    mean(v)
    error('Expected sphere with origin at (0,0,0)');
end

% coordinates of rectangle
o=v(io,:);
x=v(ix,:);
y=v(iy,:);
c=(x+y)/2; % 'center' of rectangle (not on the surface though)

ox=norm(o-x);
oy=norm(o-y);
oc=norm(o-c);

% step size along the edges of the triangle
if nx<1
    nx=ceil(1/nx)+1;
end
xstep=1/(nx-1);
if isempty(ny)
    ny=1+ceil(1/(xstep*ox/oy)); % about same step size
    ystep=1/(ny-1);
elseif ny<1
    ny=ceil(1/ny)+1;
end

maxd=max([oc,ox,oy])*1.2; % 1.2 is a safety factor to ensure we really get all nodes

% only consider nodes that are 'near'. This is a huge timesaver as later on
% we don't have to compute distances for nodes that are 'far'
d=surfing_eucldist(v',c'); % distance
nsel=find(d<=maxd);        % indices of nodes that are near
sv=v(nsel,:)';             % coordinates of nodes that are near

% allocate space for output
idxs=zeros(nx,ny);  
tic();
counter=0;
for xi=1:nx
    for yi=1:ny
        p=relpos(o,x,y,(xi-1)*xstep,(yi-1)*ystep); % position in triangle
        sidx=nearest(sv,p'); % nearest node on surface
        idxs(xi,yi)=nsel(sidx);
        
        counter=counter+1;
        if counter==0 || mod(counter,progressstep)==0 || counter==nx*ny
            surfing_timeremaining(counter/nx/ny);
        end
    end
end

function [idx,d]=nearest(xs,x) % finds the index in xs of the node nearest from x

d=surfing_eucldist(xs,x);
[dummy,idx]=min(d);
    

function y=relpos(o,x,y,a,b) % coordinates (a,b), relative to origin x and axes x and y
% TODO: fancy stuff with angles

ox=x-o;
oy=y-o;
y=o+a*ox+b*oy;



