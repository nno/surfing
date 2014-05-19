function xtfce=surfing_clustertfce(node2area,nbrs,x,dt,showprogress)
% Treshold-free cluster enhancement (TFCE) for volumes and surfaces. 
% Yields a score (e.g. t-value) to each node, that is higher if the node
%  is either part of a small but strong cluster, or a weak but large cluster
% 
% XTFCE=SURFING_CLUSTERTFCE(N2A,NBRS,X,DT)
% INPUTS:
%   N2A       Qx1 mapping from node to area for Q nodes 
%             See SURFING_SURFACEAREA
%   NBRS      QxN mapping from node to indices of neighbours
%             See SURFING_NBRS
%   X         Qx1 data to be used for TFCE; null hypothesis is that data in
%             X has zero mean and normal
%   DT        Delta of threshold (for approximation of surface integral)
%             Default is 0.1. A negative value means that only the cluster
%             sum for threshold at -DT is computed. Multiple values means 
%             that these values are taken as thresholds
% OUTPUTS:  
%   XTFCE     Qx1 TFCE values
%
% NNO Oct 2010
%
% TFCE reference: Stephen M. Smith, Thomas E. Nichols, Threshold-free
% cluster enhancement: Addressing problems of smoothing, threshold 
% dependence and localisation in cluster inference, NeuroImage, Volume 44, 
% Issue 1, 1 January 2009, Pages 83-98.
%
% See also SURFING_SURFACEAREA, SURFING_SURFACE_NBRS, SURFING_CLUSTERIZE,
% SURFING_VOLUME_NBRS

if nargin<5 || isempty(showprogress), showprogress=true; end
if nargin<4 || isempty(dt) || dt==0, dt=.1; end

if isempty(node2area)
    node2area=ones(size(nbrs,1),1);
end
nverts=numel(node2area);

if size(nbrs,1) ~= nverts, error('node2area and nbrs: different number of vertices (%d and %d)',nverts,size(nbrs,1)); end
if numel(x) ~= nverts, error('data x and nbrs: different number of vertices (%d and %d)',nverts,numel(x)); end

me=str2func(mfilename()); % for recursion

if sum(x<0) && sum(x>0)
    % separately for positive and negative values
    xtfce1=me(node2area,nbrs, x.*(x>0),dt,showprogress);
    xtfce2=me(node2area,nbrs,-x.*(x<0),dt,showprogress);
    
    % every node is either positive, negative, or zero; 
    % so it is kosher to simply take the difference
    xtfce=xtfce1-xtfce2;
    return
elseif sum(x<0) % only negative values
    xtfce=-me(node2area,nbrs,-x.*(x<0),dt,showprogress);
    return
end

% magical TFCE coeficients
B=0.5;
E=2;

% ensure that the function does not run 'almost' forever for large values
maxnsteps=1e4;

xmax=max(x(:));
if numel(dt)>1
    if ~issorted(dt) || sum(dt<=0)>0
        error('DT should be positive and sorted')
    end
    ts=dt;
elseif dt<0
    dt=-dt;
    ts=dt;
else
    nsteps=xmax/dt;
    
    if nsteps>maxnsteps;
        xmax=dt*maxnsteps;
        warning('TFCE would require more than %d steps; resetting it to %d', maxnsteps,xmax);
    end
    
    ts=dt:dt:xmax;
end

xtfce=zeros(nverts,1);
for t=ts
    clusters=surfing_clusterize(x>=t,nbrs);
    nclusters=numel(clusters);
    for k=1:nclusters
        nodeidxs=clusters{k};
        v=sum(node2area(nodeidxs))^B * t^E * dt; % TFCE formula
        xtfce(nodeidxs)=xtfce(nodeidxs)+v;
        
    end
    if showprogress
        fprintf('.');
    end
end

if showprogress
    fprintf('#\n');
end

