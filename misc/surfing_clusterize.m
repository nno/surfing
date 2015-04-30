function [c,sz,cx,xunq,xunq2cl]=surfing_clusterize(x,nbrs,n2a)
% fast depth-first clustering (supports both surface and volume)
%
% [C,SZ,CX[,XUNQ,XUNQ2CL]]=SURFING_CLUST(X,NBRS[,N2A])
% INPUTS:
%   X       Nx1 vector for N data values
%   NBRS    NxP matrix with neighbours of X (at most P per value in X),
%           zeros to indicate no neighbour.
%   N2A     Nx1 vector of size (area/volume) of each data point.
%           Default: vector with ones.
%
% OUTPUTS:
%   C       1xQ cell with cluster indices, if Q clusters are found.
%           C{k} denotes the indices in X that belong to the k-th cluster
%           Clusters are not sorted by size (but C{p}(1) < C{q}(1) if p<q)
%   SZ      Qx1 vector with size of each cluster
%   CX      Qx1 vector with CX{k} the value of X in the k-th cluster
%   XUNQ    Sx1 vector with unique elements in X
%   XUNQ2CL 1xS cell with cluster indices of unique elements.
%           Thus, CX(XUNQ2CL{k}) is a vector with only the value XUNQ(k)
%
% Two values X(j) and X(k) are in the same cluster if X(j)==X(k), and if
% they are connected (where connectedness is the transitive closure of
% neighbourness, which in turn is assumed to be symmetric).
% Values of zero and NaN in X are ignored for clusters.
%
% Example: suppose x is a map with z-scores, and clusters are to be found
% for which the z score is significant at |z|>=2. It is assumed that nbrs is
% properly defined. Define xsign=(X>=2)-(X<=-2); so that Xsign has the
% value +/-1 for signifcantly above/below chance z scores. Then
%
%   [c,sz,cx,xunq,xunq2cl]=surfing_clust(xsign,nbrs)
%
% finds the required clusters required. Moreover, xunq2cl{1} and xunq2cl{2}
% contains the cluster indices significantly below and above chance.
%
% See also: surfing_surfacearea, surfing_surface_nbrs, surfing_volume_nbrs
%
% NNO Sep 2010, updated Jan 2011

[n,p]=size(nbrs);
[n2,p2]=size(x);

if n2~=n || p2 ~=1
     error('illegal input size');
end

if nargin<3
    n2a=ones(n,1);
elseif numel(n2a) ~= n
    error('n2a has illegal size');
end

q=zeros(n,1); % queue of nodes to add to current cluster
qend=0;
qnext=1; % first free position

clstart=zeros(n+1,1); % start index of i-th cluster. One bigger for ease in setting output
clcount=0; % total number of clusters found so far
clidxs=zeros(n,1); % label for each cluster (positive integers)
inqueue=false(n,1); % indicates if an element is in the queue

xpos=1; % first non-visited node

% go over all nodes from 1 to n. xpos is the candidate node for a new
% cluster, and we add its neighbours in a breadth-first search
while xpos<=n
    if clidxs(xpos)>0 % is already part of a cluster, continue
        xpos=xpos+1;
        continue;
    end

    xval=x(xpos); % x value of candidate cluster

    if xval~=0 && ~isnan(xval)

        % found an element for a new cluster

        q(qnext)=xpos; % first element in the queue
        qend=qnext; % last element of the queue
        clcount=clcount+1; % number of clusters
        clstart(clcount)=qnext; % first index of this cluster
        clidx=clcount; % index for current clustter

        inqueue(xpos)=true;

        % now visit all neighbours of q(qnext), iteratively, until the queue is
        % exhausted
        while true
            qxpos=q(qnext); % position in queue
            clidxs(qxpos)=clidx; % assign cluster index

            for j=1:p % check all nbrs
                xnbr=nbrs(qxpos,j);
                % neighbour should have positive index, not queued yet, and
                % have the same value as the first element of the current
                % cluster
                if xnbr>0 && ~inqueue(xnbr) && x(xnbr)==xval % neighbour and correct value
                    qend=qend+1;
                    q(qend)=xnbr;
                    inqueue(xnbr)=true; % mark as cluster that is to be visited, not again added to the queue
                end
            end

            % get ready for next element in queue
            qnext=qnext+1;

            % if queue is exhausted, break out and try to find a new cluster
            if qnext>qend
                break;
            end
        end
    end
    xpos=xpos+1;
end

% allocate space for output
c=cell(1,clcount); % cluster index
sz=zeros(1,clcount); % cluster size
cx=zeros(1,clcount); % cluster value (from x)

% set extra boundary for last cluster
clstart(clcount+1)=qend+1;

% assign cluster index, size, and value
for k=1:clcount
    qk=q(clstart(k):(clstart(k+1)-1));
    c{k}=qk;
    sz(k)=sum(n2a(qk));
    cx(k)=x(qk(1));
end

if nargout>4
    % this could be time-consuming if there are many different values in x
    [xunq,is]=unique(cx);
    nis=numel(is); % number of unique elements
    xunq2cl=cell(1,nis);
    for k=1:nis
        xunq2cl{k}=find(cx==cx(is(k)));
    end
end






