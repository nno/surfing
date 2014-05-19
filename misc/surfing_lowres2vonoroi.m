function [vonl2h,vonh2l,mindist]=surfing_lowres2vonoroi(c1,c2,f2,r)
% constructs a Vonoroi mapping on a high resolution surface, using a low
% resolution surface as seed points and geodesic distance measure. 
%
% [L2H,H2L,D]=SURFING_LOWRES2VONOROI(CLOW,CHIGH,FHIGH[,R]) 
% INPUTS:
%   CLOW:   3xNLOW  coordinates for N nodes from low  resolution surface
%   CHIGH:  3xNHIGH "                          " high "                "
%   FHIGH   3xP vertex indices for P faces (base1)
% OPTIONAL INPUT:
%   R:      Initial radius on high resolution surface; default is R=6.
%           Because the necessary Euclidian distance for search regions 
%           cannot be computed beforehand, a certain radius R is chosen 
%           initially and increased until for all nodes the minimum 
%           distance was found.
%           A negative value of R will use radius -R, but not increase R,
%           meaning that certain values in D may be Inf.
% OUTPUTS:
%   L2H:    NLOWxQ node mapping from low to high resolution surface nodes.
%           At most Q ndoes in CHIGH are associated with a node in CLOW
%   H2L:    NHIGHx1 node mapping from high to low resolution surface nodes.
%   D:      NHIGHx1 geodesic distances from a node in CHIGH and the nearest 
%           corresponding node in CLOW. 
%
% It is assumed that the low and high resolution surfaces were constructed
% using AFNI's or Freesurfer's mapicosehedron (see SURFING_MAPLOW2HIRES)
%
% NNO May 2010
%
% See also SURFING_MAPLOW2HIRES, SURFING_NODEIDXS2FACEIDXS,
% SURFING_CIRCLEROI


if nargin<4
    r=6; 
end
rinc=6;  %increment of radius
rmult=1; %multiply radius

returnimmediately=r<0;
if returnimmediately
    r=-r;
end

l2h=surfing_maplow2hires(c1,c2);     %find mapping from low to hi res surface
n2f2=surfing_nodeidxs2faceidxs(f2);  % find mapping from nodes to faces on hi res surface

nverts1=size(c1,2);
nverts2=size(c2,2);

while true
    tic();
    fprintf('surfing_lowres2vonoroi: using radius=%.1f\n',r);
    mindist=repmat(Inf,nverts2,1);
    vonh2l=zeros(nverts2,1);

    for k=1:nverts1
        ci=l2h(k);
        [idxs,fff,dist]=surfing_circleROI(c2,f2,ci,r,n2f2);

        msk=dist<=mindist(idxs);
        vonh2l(idxs(msk))=k;
        mindist(idxs(msk))=dist(msk);
    end
    
    infcount=sum(isinf(mindist));
    toc();
    if infcount==0 || returnimmediately
        fprintf('Completed using radius %.1f\n',r);
        break;
    else
        rnew=r*rmult+rinc;
        fprintf('Did not complete %d / %d nodes (%.1f%%) using radius %.1f, incrementing radius to %.1f \n',infcount,nverts2,100*infcount/nverts2,r,rnew);
        r=rnew;
    end
end

vonl2h=surfing_invertmapping(vonh2l);
