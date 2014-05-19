function [sc, sf, si, nsel, fsel] = surfing_subsurface(c, f, i, r, n2f)
% Generates a sub-surface with nodes from a surface that are within a 
% certain radius (Euclidian distance) from a center node.
%
% [SV, SF, SI, NSEL, FSEL] = SURFING_SUBSURFACE(C, F, I, R[, N2V]) 
% INPUTS:
%   C:         3xN coordinates for N nodes
%   F:         3xP vertex indices for P faces, base1
%   I:         single integer in range 1:N indicating the center node index
%   R:         radius used to select nodes. Nodes within an Euclidian
%              distance from the center node less than R are returned.
% OPTIONAL INPUT:
%   N2F:       NxQ mapping from node indices to face indices, if each node 
%              is contained in at most Q faces (typically Q==6). If omitted, 
%              this mapping is generated using SURFING_NODEIDXS2FACEIDXS.
% OUTPUTS:
%   SC,SF      coordinates and faces of the subsurface. 
%   SI:        index of center node relative to subsurface; C(I,:)==SC(SI,:)
%   NSEL,FSEL: indices of nodes and faces selected from C and F.
%
% Experimental features: if R is empty and I contains multiple indices, then a
% subsurface is returned containing at least vertices I and the faces
% surrounding I.
%
% Note that this function may return some nodes that are actually further
% away from node I than radius R. The reason is that this function returns
% a sub-surface that contains all nodes within distance R, then selects all
% the faces that contain at least one of these nodes, and finally returns a
% sub-surface with all the nodes that are contained in at least one of
% these faces.
%
% NNO May 2009, May 2010
%
% See also SURFING_CIRCLEROI, SURFING_NODEIDXS2FACEIDXS


% check input size, transpose if necessary
trn=size(c,1) ~= 3; if trn, c=c'; end
trf=size(f,1) ~= 3; if trf, f=f'; end

if size(c,1) ~= 3 || size(f,1) ~= 3, error('Coordinates and faces should be 3xP and 3xQ'); end

% see if we have to compute node to face mapping on the fly
showwarning=nargin<5 || isempty(n2f); % always show a warning if n2f is not specified
if showwarning,n2f=[];end
if isempty(n2f);
    if showwarning
        warning('SURFING:notoptimal',sprintf('n2f (node to voxel indices mapping) not specified.\nWill compute on the fly, but this can be slow.\nConsider running SURFING_NODEIDXS2VOXELIDXS first'));
    end
    n2f=surfing_nodeidxs2faceidxs(f);
end

if size(c,2) ~= size(n2f,1), error('size of N2F and V do not match'); end

if isempty(r)
    vmsk=i(:); % take the node indices directly [EXPERIMENTAL]
elseif numel(i)==1
    % compute euclidian distance from center node
    vmsk=surfing_eucldist(c(:,i),c)<=r;
else
    error('illegal input for i and r')
end

fidxs=n2f(vmsk,:); %indices of faces corresponding to node mask
funq=unique(fidxs(fidxs>0)); %unique face indices

fsel=f(:,funq);             % nodes corresponding to face indices
[nsel,p,q]=unique(fsel(:)); % find unique nodes
sc=c(:,nsel);                    % coordinates of selected nodes
sf=reshape(q,3,numel(funq));     % new topology with limited number of nodes
si=find(nsel==i); %can this be done more efficiently?

% transpose back if necessary
if trn, sc=sc'; end
if trf, sf=sf'; end