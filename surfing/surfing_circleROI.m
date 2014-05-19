function [coordidx,D,scoords,vORr]= surfing_circleROI(coords,faces,centervertex_idx,radius,distancemetric,n2f) 
% Creates a circular region of interest (ROI) on the surface given a radius
%
% [coordidx,D,scoords]=surfing_circleROI(coords,faces,centervertex_idx,radius,distmetric,n2f)
% INPUTS:
%   coords:     3xN coordinates for the N vertices
%   faces:      3xP vertex indices for the P triangualar faces (1-based indices)
%   centeridx:  the index of the center vertex 
%   radius:     the radius of the circle; use [R C] to select C nodes with
%               initial radius R
% OPTIONAL INPUTS
%   distmetric: distancemetric: 'euclidian' or 'dijkstra' or 'geodesic' (default)
%   n2f:        for faster computation, a to face mapping N2V (NxM, if each node is contained 
%               in M faces at most); see SURFING_NODEIDXS2FACEIDXS. 
%               This option is required only for 'dijkstra' and
%               'geodesic_legacy'.
% OUTPUT: 
%   coordidx:   1xK vector with the K vertex indices that are within distance RADIUS from the 
%               center vertex
%   D:          1xK vector of the distances from center vertex
%   scoords:    3xK matrix of coordinates of the selected vertices
%
% If distmetric is omitted, then a geodesic distance argument is assumed; N2F can 
% still be passed as the sixth argument.
% 
% Computation of geodesic distances uses the Fast Marching toolbox by 
% Gabriel Peyre (2008), http://www.ceremade.dauphine.fr/~peyre/download/
% 
% NNO,TW,JD May 2010
%
% See also SURFING_NODEIDXS2FACEIDXS, SURFING_SUBSURFACE,
% PERFORM_FAST_MARCHING_MESH

if nargin<5
    n2f=[];
	distancemetric='geodesic';
elseif isnumeric(distancemetric)
    n2f=distancemetric;
    distancemetric='geodesic';
end

% NNO Apr 2011 support for fixed number of nodes (variable radius)
variableradius=numel(radius)>1;
if variableradius
    radiusgrow=1.5; 
    radiuswarn=200;
    if radius(1)==0
        radius(1)=10;
    end
end


% check size of coords and faces, transpose if necessary
if size(coords,1)~=3, coords=coords'; end
if size(faces,1)~=3 , faces=faces'; end
if ~isequal([size(coords,1),size(faces,1)],[3 3])
    error('Expected three dimensional coords and faces matrices');
end

while true
    switch distancemetric
        case 'euclidian'
            D=surfing_eucldist(coords(:,centervertex_idx),coords)';  
            vidxs=(1:size(coords,2))'; % all vertex indices

        case 'geodesic'
            % as of June 2013 this is obsolete. Releases before this date
            % used the code below
            [sv, sf, si,vidxs, fidxs]=surfing_subsurface(coords, faces, centervertex_idx, radius(1), n2f); % construct correct sub surface

            % this requires the Fast Marching toolbox (Peyre)
            [D,S,Q] = perform_fast_marching_mesh(sv, sf, si);     
            
        case 'dijkstra'
            [sv, sf, si,vidxs, fidxs]=surfing_subsurface(coords, faces, centervertex_idx, radius(1), n2f); % construct correct sub surface
            D=surfing_dijkstradist(sv,sf,si,radius(1));
            
        otherwise
            error('Unknown distance metric %s, use ''geodesic'' or ''dijkstra'' or ''euclidian''',distancemetric);
    end
    
    nodemask=D<=radius(1); % selected nodes

    if ~variableradius || sum(nodemask)>=radius(2)  % big enough      
        break;
    end
    
    % too small; increase radius and try again
    radius(1)=radius(1)*radiusgrow;
    if radius(1)>radiuswarn
        warning('Radius has become really big: %d mm', radius(1));
    end
end

% set final resutls
if variableradius
    [dummy,is]=sort(D,'ascend');
    issel=is(1:radius(2));
    D=D(issel)';
    selidxs=issel;
    vORr=D(end);
else
    selidxs=find(nodemask);
    D=D(selidxs)';
    vORr=numel(selidxs);
end
coordidx=vidxs(selidxs)';
scoords=coords(:,coordidx);

return
