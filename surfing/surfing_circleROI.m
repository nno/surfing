function [coordidx,D,scoords,vORr]= surfing_circleROI(coords,faces,centervertex_idx,radius,distancemetric,n2f) 
% Creates a circular region of interest (ROI) on the surface given a radius
%
% [coordidx,D,scoords]=surfing_circleROI(coords,faces,centervertex_idx,...
%                                           radius,distmetric,n2f)
% INPUTS:
%   coords:     3xN coordinates for the N vertices
%   faces:      3xP vertex indices for the P triangualar faces 
%               (1-based indices)
%   centeridx:  the index of the center vertex 
%   radius:     the radius of the circle:
%               - use R to select all nodes within radius R (typically in mm)
%               - use [R C] to select C nodes (with initial radius R)
%               - use [R Inf A] to select the nearest nodes nodes whose 
%                 area together is less than or equal to area A 
%                 (typically in mm^2)
% OPTIONAL INPUTS
%   distmetric: distance metricL 'euclidean', 'dijkstra', 
%               or 'geodesic'  (default)
%   n2f:        for faster computation, a to face mapping N2V (NxM, 
%               if each node is contained in M faces at most); see 
%               SURFING_NODEIDXS2FACEIDXS. 
%               If omitted this mapping is computed on the fly, which is
%               more time consuming if this function is called multiple
%               times (use case: searchlight analysis)
% OUTPUT: 
%   coordidx:   1xK vector with the K vertex indices that are within 
%               distance RADIUS from the center vertex
%   D:          1xK vector of the distances from center vertex
%   scoords:    3xK matrix of coordinates of the selected vertices
%   vORr        if radius is of the form
%               R         : the number of vertices selected
%               [R C]     : the final radius
%               [R Inf A] : the final surface area
%
% If distmetric is omitted, then a geodesic distance argument is assumed; 
% N2F can still be passed as the sixth argument.
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
    radiusmax=200;
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

% node areas
based_on_area=numel(radius)==3;
if based_on_area
    node2area=surfing_surfacearea(coords,faces,n2f);
end

skip_node=any(isnan(coords(:,centervertex_idx)));
if skip_node
    coordidx=zeros(1,0);
    D=zeros(1,0);
    scoords=zeros(3,0);
    return
end

while true
    switch distancemetric
        % earlier versions had a typo in euclidean, so support
        % both the name with the typo and without
        case {'euclidian','euclidean'}
            D=surfing_eucldist(coords(:,centervertex_idx),coords)';  
            vidxs=(1:size(coords,2))'; % all vertex indices

        case 'geodesic'
            % construct correct sub surface
            [sv, sf, si,vidxs, fidxs]=surfing_subsurface(coords, faces, ...
                                        centervertex_idx, radius(1), n2f); 

            % this requires the Fast Marching toolbox (Peyre)
            [D,S,Q] = perform_fast_marching_mesh(sv, sf, si);     
            
        case 'dijkstra'
            % construct correct sub surface
            [sv, sf, si,vidxs, fidxs]=surfing_subsurface(coords, faces, ...
                                        centervertex_idx, radius(1), n2f); 
            D=surfing_dijkstradist(sv,sf,si,radius(1));
            
        otherwise
            error(['Unknown distance metric %s, use ''geodesic'', '...
                    '''dijkstra'' or ''euclidean'''],distancemetric);
    end
    
    nodemask=isfinite(D) & D<=radius(1); % selected nodes

    if ~variableradius || sum(nodemask)>=radius(2)  % big enough      
        break;
    end
    
    % check area
    if based_on_area && sum(node2area(vidxs(nodemask))) >= radius(3)
        break;
    end
    
    % too small; increase radius and try again
    radius(1)=radius(1)*radiusgrow;
    if radius(1)>=radiusmax
        if isinf(radius(1))
            error(['unable to select %d nodes around node %d, maximum '...
                        'count is %d'], ...
                        radius(2) , centervertex_idx, sum(nodemask));
        else
            radius(1)=Inf;
            % and try on next iteration try to get all nodes
        end
    end
end

% set final resutls
if based_on_area
    [dummy,is]=sort(D,'ascend');
    cumarea=cumsum(node2area(vidxs(is)));
    msk=cumarea<=radius(3);
    selidxs=is(msk);
    D=D(selidxs);
    vORr=max(cumarea(msk));
elseif variableradius
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
