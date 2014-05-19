function [qr_pth, fi, fj, cp]=surfing_node_half_collapse(v,f,vq,vri,pths,v2f,f_cross)
% determine surface changes for a half-collapse of a surface node
%
% [qr_pth, fi, fj, cp]=surfing_node_half_collapse(v,f,vq,vri,pths,v2f,f_cross)
%
% Inputs:
%   v           Px3 coordinates for P nodes
%   f           Qx3 node indices for Q faces
%   vq          Coordinate index of pivot node (that is removed)
%   vri         Path index of node into which is collapsed, assuming that
%               pths(vq,vri) is non-zero
%   pths        PxM path indices for each node in v, if each node as a path
%               of at most M nodes (see surfing_surface_simple_nodes)
%   v2f         Mapping from vertices to faces (see 
%               surfing_nodeidxs2faceidxs)
%   f_cross     Qx3 cross product of two vertices of each face, i.e. 
%               if f(k,:)==[a,b,c] then it is the cross product of
%               vertices a-b and a-c.
% Returns:
%   qr_pth      Indices around node vq and pths(vq,vri)
%   fi          1x2 face indices to be removed 
%   fj          1xR face indices to be kept
%   cp          Rx1 Surface ratios (original versus new) of faces in fj
%   
% Notes:
%   - this function returns empty output if:
%     * node vq is not a simple node
%     * pths(vq,vri)==0
%     * there are not exactly two faces that contain both nodes
%       vq and pths(vq,vri)
%     * collapsing node vq into node pths(vq,vri) would flip a triangle
%   - to subsample a surface, consider using surfing_subsample_surface
%
% Example:
%  - consider the following surface, with arabic numerals indicating node
%    indices and roman numerals face indices:
%
%          2 ----- 7
%         / \     / \
%        / I \ VI/V  \
%       /     \ /     \
%      3 ----- 1 ----- 6
%       \     / \     /
%        \ II/III\ IV/
%         \ /     \ /
%          4 ------5
%
%    this surface can be described by
%    >> v=[0 -1 -2 -1 1 2 1;0 -2 0 2 2 0 -2;0 0 0 0 0 0 0]';
%    >> f=[1 1 1 1 1 1; 2 3 4 5 6 7; 3 4 5 6 7 2]';
%
%    and the other properties can be computed by
%    
%    >> pths=surfing_surface_simple_nodes(v,f);
%    >> v2f=surfing_nodeidxs2faceidxs(f');
%    >> a=v(f(:,1),:)';
%    >> b=v(f(:,2),:)';
%    >> c=v(f(:,3),:)';
%    >> f_cross=cross(a-b,a-c)';
%
%    Only node 1 can be collapsed because it is the only simple node.
%    Collapsing node 1 into node 6 would yield the following (with all
%    lines being straight, which is not shown properly due to ASCII
%    limitations):
%
%          2 ----- 7
%         /  \  VI  \
%        / I   ----- \
%       /            \\
%      3 ------------- 6
%       \            //
%        \ II  ----- /
%         \  / III  /
%          4 ------5
%
%    >> vq=1;
%    >> vr=6
%    >> vri=find(pths(vq,:)==vr)
%    >> [qr_pth, fi, fj, cp]=surfing_node_half_collapse(v,f,vq,vri,pths,v2f,f_cross)
% 
%    qr_pth =
% 
%      2     3     4     5     6     7     1     6
% 
% 
%    fi =
% 
%      4     5
% 
% 
%    fj =
% 
%      1     2     3     6
% 
% 
%    cp =
% 
%      2
%      2
%      1
%      1
%    
%     Thus qr_pth contains the node indices around node 1 and node 6 (all
%     nodes in this case), if contains the indices of faces containing node
%     1 and node 4 (faces IV and V), fj contains the indices of faces
%     that contain node 1 but node node 4 (faces I, II, III and VI), and cp 
%     is ratio between surfaces of the kept faces if the node collapse 
%     would be performed (faces I and II double in size, faces III and VI
%     maintain their size)
%
%
% See also: surfing_surface_simple_nodes, surfing_nodeidxs2faceidxs,
%           surfing_subsample_surface
% 
% NNO May 2014



fi=[];
fj=[];
qr_pth=[];
cp=[];

while true
    pth=pths(vq,:);
    pth=pth(pth>0);
    
    n=numel(pth);

    if n==0 || vri>n
        break;
    end
    
    vr=pth(vri);
           
    fq=v2f(vq,:);
    fr=v2f(vr,:);
    fq=fq(fq>0);
    fr=fr(fr>0);

    m=bsxfun(@eq,fq',fr);
    m_q=sum(m,2)>0;

    if all(m_q)
        break
    end

    if sum(m_q)~=2
        break
    end
    assert(sum(m_q)==2)

    f_q_orig=f(fq(~m_q),:);
    cross_orig=f_cross(fq(~m_q),:);

    f_q_new=f_q_orig;
    f_q_new(f_q_new==vq)=vr;

    a=v(f_q_new(:,1),:)';
    b=v(f_q_new(:,2),:)';
    c=v(f_q_new(:,3),:)';
    cross_new=cross(a-b,a-c)';
    
    inner=sum(cross_new.*cross_orig,2);
    if ~all(inner>0)
        break;
    end
    cp=sqrt(sum(cross_new.^2,2)./sum(cross_orig.^2,2));
    
    qr_pth=[pth pths(vr,:) vq vr];
    qr_pth=qr_pth(qr_pth>0);
    
    fi=fq(m_q);
    fj=fq(~m_q);
    break;
end
    
