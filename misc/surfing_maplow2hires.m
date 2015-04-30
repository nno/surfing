function [low2hi,ds]=surfing_maplow2hires(clow, chigh, epsilon)
% finds mapping between nodes from low and high res surfaces
%
% [L2H,D]=SURFING_MAPLOW2HIRES(CLOW,CHIGH[,EPSILON])
% INPUTS:
%   CLOW:     3xNLOW  coordinates for N nodes from low  resolution surface
%   CHIGH:    3xNHIGH "                          " high "                "
%   EPSILON:  Scalar value of maximum difference between nodes in CLOW and
%             CHIGH (default: 1e-4)
% OUTPUTS:
%   L2H:    NLOWx1 node mapping from low to high resolution surface nodes.
%           If L2H(I)==J, then node J on CHIGH is the nearest node to node
%           J on CLOW.
%   D:      NLOWx1 eulcidian distance from node I in CLOW and node D(I) in
%           CHIGH.
%
%
% NNO Aug 2009, updated May 2010.

if nargin<3
    epsilon=1e-4;
end


% number of vertices
if size(clow,1) ~= 3 || size(chigh,1) ~= 3
    error('Coordinates should be 3xN');
end

nlow=size(clow,2);
nhi=size(chigh,2);

if nlow>nhi
    warning('Lo res surface has more (%d > %d) nodes than hi res',...
                nlow,nhi);
elseif nlow==nhi && isequal(clow,chigh)
    % little optimization: same surfaces, return identity mapping
    low2hi=(1:nhi)';
    ds=zeros(nhi,1);
    return
end

low2hi=zeros(nlow,1);
ds=zeros(nlow,1);

clock_start=clock();
prev_msg='';
progress_step=100;

% progress message
msg=sprintf('%d -> %d nodes', nhi, nlow);

if false
    step_size=100;
    pos=0;
    for k=1:step_size:nlow
        idxs=pos+(1:step_size);

        hids=surfing_eucldist(clow(:,idxs),chigh);

        % find nearest high res vertex
        [d,idx]=min(hids,[],2);

        if any(d>epsilon)
            i=idxs(idx);
            error('Low res node %d at %.3f>%.3f from nearest node',...
                        i,d(i),epsilon);
        end

        low2hi(idxs)=idx;
        ds(idxs)=d;

        pos=pos+step_size;
        if pos+step_size>nlow
            step_size=nlow-pos;
        end

        prev_msg=surfing_timeremaining(clock_start,k/nlow,msg,prev_msg);
    end


else


    for k=1:nlow
        % distance from low res vertex k to all high res vertices
        hids=surfing_eucldist(clow(:,k),chigh);

        % find nearest high res vertex
        [d,idx]=min(hids);


        if d>epsilon
            error('Low res node %d has smallest distance %.3f to high res node %d, exceeding epsilon=%.3f',k,d,idx,epsilon);
        end

        low2hi(k)=idx;
        ds(k)=d;

        if mod(k,progress_step)==0
            prev_msg=surfing_timeremaining(clock_start,k/nlow,msg,prev_msg);
        end
    end
end
