function [vv,ff]=surfing_subsample_surface(v,f,niter,min_ratio,...
                                                        progress_step)
% removes a subset of nodes to decrease the node density of a surface
%
% [vv,ff]=surfing_subsample_surface(v,f[,niter[,min_ratio]])
%
% Inputs:
%   v           Px3 coordinates for P nodes
%   f           Qx3 node indices for Q faces
%   niter       Number of iterations (default: 1)
%   min_ratio   Minimal surface ratio between original and replacement
%               faces; default .2
%
% Outputs:
%   vv          PPx3 coordinates of PP nodes (PP<=P)
%   ff          QQx3 node indices for QQ faces (QQ<=Q)
%
% Example:
%  - consider the following surface, with arabic numerals indicating node
%    indices and roman numerals face indices:
%
%          2 ----- 7 ----- 10
%         / \     / \     / \
%        / I \ VI/V  \ X / IX\
%       /     \ /     \ /     \
%      3 ----- 1 ----- 6 ----- 9
%       \     / \     / \ VIII/
%        \ II/III\ IV/VII\   /
%         \ /     \ /     \ /
%          4 ----- 5 ----- 8
%
%    this surface can be described by
%    >> v=[0 -1 -2 -1 1 2 1 3 4 3;0 -2 0 2 2 0 -2 2 0 -2;zeros(10,1)]';
%    >> f=[1 1 1 1 1 1 5 8 6 6;2 3 4 5 6 7 8 9 9 10;3 4 5 6 7 2 6 6 10 7]';
%
%    Subsampling this surface for two iterations gives:
%    >> [vv,ff]=surfing_subsample_surface(v,f,2)
%
%    vv =
%
%     -1    -2     0
%     -2     0     0
%     -1     2     0
%      1     2     0
%      1    -2     0
%      3     2     0
%      4     0     0
%      3    -2     0
%
%
%    ff =
%
%      1     2     3
%      1     3     4
%      1     4     5
%      4     8     5
%      6     7     4
%      4     7     8
%
%   corresponding to (with also edge 4-7 being a straight line)
%
%           1 ----- 2 ----- 8
%          /|\      |      / \
%         / | \ III | IV  /   \
%        /  |  \    |    /     \
%       2 I |   \   |   / VI    7
%        \  |    \  |  /      //
%         \ | II  \ | / ----- /
%          \|      \|// V    /
%           3 ----- 4 ----- 5
%
%   which shows a reduction by two nodes and 4 faces.
%
% Notes:
%  - this function uses a repeated application of
%    node half-collapse; for details see surfing_node_half_collapse
%  - only 'simple' nodes are replcated; for details see
%    surfing_surface_simple_nodes
%
% See also: surfing_surface_simple_nodes, surfing_node_half_collapse
%
% NNO May 2014

if nargin<5 || isempty(progress_step)
    progress_step=5000;
end

if nargin<4 || isempty(min_ratio)
    min_ratio=.2;
end

if nargin<3 || isempty(niter)
    niter=1;
end

nv=size(v,1);

if niter==0
    vv=v;
    ff=f;
    return
elseif niter>1
    for iter=1:niter

        [vv,ff]=surfing_subsample_surface(v,f,1,min_ratio,progress_step);
        nv=size(v,1);
        nvv=size(vv,1);
        if nvv==nv
            break;
        end
        v=vv;
        f=ff;
    end
    return
end

nf=size(f,1);

% get 'simple' nodes - only those are potentially removed
[pths,is_simple]=surfing_surface_simple_nodes(v,f,progress_step);
npths=sum(pths>0,2);
assert(~any(npths>0 & npths<3));

% compute cross product for two vertices of each triangle
v2f=surfing_nodeidxs2faceidxs(f');
a=v(f(:,1),:)';
b=v(f(:,2),:)';
c=v(f(:,3),:)';
f_cross=cross(a-b,a-c)';

% keep track of which vertices and vertices are to be kept
f_keep=true(nf,1);
v_keep=true(nv,1);

% mapping from old to new node assignment
old2new=zeros(nv,1);

% keep track of which nodes were visited
visited=false(nv,1);

% a counter to pseudo-randomly remove nodes
counter=0;

% define a queue of nodes to be considered
queue=zeros(nv+1,1);
queue_start=0;
queue_end=0;

% show progress to the user
clock_start=clock();
prev_msg='';
nv_removed=0;

while ~all(visited | ~is_simple)
    for npth=3:max(npths)

        % get candidates with path around it of 'col' nodes
        candidates=find(npths==npth & ~visited);
        if isempty(candidates)
            continue;
        end

        % add candidates to queue
        ncandidates=numel(candidates);
        queue(queue_start+(1:ncandidates))=candidates;
        queue_end=queue_end+ncandidates;

        % loop over node candidates
        while queue_start<queue_end
            % pop a new node candidate from the queue
            queue_start=queue_start+1;
            vq=queue(queue_start);

            % show progress
            if progress_step && ...
                    (queue_start==1 || mod(queue_start, progress_step)==0)
                msg=sprintf('%d nodes: %.1f%% removed', ...
                                        nv, nv_removed*100/nv);
                prev_msg=surfing_timeremaining(clock_start,...
                                            queue_start/nv,msg,prev_msg);
            end

            % if already visited then ignore this node
            if visited(vq)
                continue;
            end

            % mark node as visited
            visited(vq)=true;

            for j=1:npth
                % select pseudo-randomly a node 'vr' neighboring node 'vq'
                counter=counter+1;
                vri=mod(counter-1,npth)+1;
                vr=pths(vq,vri);
                if visited(vr)
                    continue;
                end

                % faces containing nodes vq and vr
                fq=v2f(vq,:);
                fr=v2f(vr,:);

                % get rid of empty elements
                fq=fq(fq>0);
                fr=fr(fr>0);

                % see which faces are common in both nodes
                m=bsxfun(@eq,fq',fr);

                % mask for vq that is 'true' wherever a face in fq is also
                % contained in fr
                m_q=sum(m,2)>0;

                if all(m_q)
                    continue
                end

                % faces in node vq not in common with node vr; these
                % faces are kept
                q_f_keep=fq(~m_q);

                % get cross product for each of the faces
                cross_orig=f_cross(q_f_keep,:);

                % replace node vq by node vr in these faces
                f_q_new=f(q_f_keep,:);
                f_q_new(f_q_new==vq)=vr;

                % compute cross product after replacement
                a=v(f_q_new(:,1),:)';
                b=v(f_q_new(:,2),:)';
                c=v(f_q_new(:,3),:)';
                cross_new=cross(a-b,a-c)';

                % compute inner product between original and new cross
                % products; negative values indicate an (undesirable)
                % triangle flip
                inner=sum(cross_new.*cross_orig,2);

                % if any replacement caused a triangle flip then
                % do not continue
                if ~all(inner>0);
                    continue
                end

                % compute ratios between original and new triangle areas.
                % if any ratio is too small then do not continue
                ratios=sqrt(sum(cross_new.^2,2)./sum(cross_orig.^2,2));
                if any(ratios<min_ratio)
                    continue
                end

                % all good, apply the changes to the faces
                f(q_f_keep,:)=f_q_new;

                % store mapping from old to new
                old2new(vq)=vr;

                % set node and faces to be removed
                v_keep(vq)=false;
                f_keep(fq(m_q))=false;

                % keep track of how many nodes were removed
                nv_removed=nv_removed+1;

                % set neighbors of nodes vq and vr as visited, so that
                % they are not going to be removed.
                qr_nbrs=[pths(vq,:) pths(vr,:)];
                visited(qr_nbrs(qr_nbrs>0))=true;
            end
        end
    end
end

% remove nodes and faces marked for removal
vv=v(v_keep,:);
[unq,p,q]=unique(f(f_keep,:));
ff=reshape(q,[],3);

% update progress
if progress_step
    msg=sprintf('%d nodes: %.1f%% removed', nv, nv_removed*100/nv);
    surfing_timeremaining(clock_start,1,msg,prev_msg);
end

% remove duplicate triangles
[sff,i]=sortrows(sort(ff')');
m=all(diff(sff)==0,2);
ff=ff(i([true; ~m]),:);

% ensure all is kosher; if not give a warning message
msg=surfing_check_surface(vv,ff,false,2);
if ~isempty(msg)
    warning(msg);
end
