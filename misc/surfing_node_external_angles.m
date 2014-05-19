function [theta_ext, nbrs]=surfing_node_external_angles(v,f)
% compute external angles of neighbors for each node
%
% [theta_ext, nbrs]=surfing_node_external_angles(v,f)
% 
% Inputs:
%   v              Px3 node coordinates
%   f              Qx3 face indices
%
% Outputs:
%   theta_ext      PxM angle (in radians), if each node has at most M
%                  neighboring nodes
%   nbrs           PxM indices of neighboring nodes
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
%    This surface has 7 nodes and 6 faces. If a-b-c denotes the angle at 
%    node b made by the edges a-b and b-c, then node 1 has 6 external 
%    angles (as it has 6 neighbors), namely 2-3-4, 3-4-5, 5-6-7, 6-7-2, 
%    and 7-2-3; all other nodes have 3 external angles.
%    * nbrs(i,:) contains the indices of nodes that are a neighbor of
%      node i; each entry for i>1 has three values (three neighbors) and 
%      zeros elsewhere.
%    * theta contains the corresponding angles. For example the neighbors
%      of node 6 (6-th row of nbrs) are 5, 1, and 7. The corresponding 
%      external angle 4-5-6 is ~0.92 radians, 5-1-7 is ~2.03 radians, and
%      6-7-1 is ~0.92 radians.
%
%    Using this function gives
%    >> [theta_ext,nbrs]=surfing_node_external_angles(v,f)
%
%    theta_ext =
% 
%     2.2143    2.0344    2.0344    2.2143    2.0344    2.0344
%     1.1071    2.0344    1.1071         0         0         0
%     0.9273    2.2143    0.9273         0         0         0
%     1.1071    2.0344    1.1071         0         0         0
%     1.1071    2.0344    1.1071         0         0         0
%     0.9273    2.2143    0.9273         0         0         0
%     1.1071    2.0344    1.1071         0         0         0
% 
% 
%    nbrs =
% 
%      3     4     5     6     7     2
%      3     1     7     0     0     0
%      2     1     4     0     0     0
%      3     1     5     0     0     0
%      4     1     6     0     0     0
%      5     1     7     0     0     0
%      6     1     2     0     0     0
%
% Notes:
%   - To convert angles to degrees, multiply by (180/pi).
%
% NNO Apr 2014

[theta_int,unused,v2f]=surfing_node_internal_angles(v,f);
nv=size(v,1);
nf=size(f,1);

nbrs=surfing_surface_nbrs(f);
n_nbrs_max=size(nbrs,2);
n_f_max=size(v2f,2);
theta_ext=zeros(nv,n_nbrs_max);

for col=1:n_nbrs_max
    j=nbrs(:,col); % neighbor of each node
    mj=find(j>0);        % node mask that have a neighbor
    
    fmi=v2f(mj,:);    % faces of source node            } masked 
    fmj=v2f(j(mj),:); % faces of each neighbor node j   } by mj
    
    % for each column in fmj, see for which rows there is a matching
    % element in the corresponding row of fmi. In these rows the faces
    % have a node in common. 
    for col_=1:n_f_max
        % consider each face that contains neighbor node j
        
        % nodes (masked by mj) for which the neighbor j has a face in
        % common in the 'col_'-th column
        %mij=any(bsxfun(@eq,fmj(:,col_),v2f(mj,:)),2);
        m=any(bsxfun(@eq,fmj(:,col_),fmi),2);
        
        mm=mj(m);
        
        theta_ext(mm,col)=theta_ext(mm,col)+abs(theta_int(j(mm),col_));
    end
    
        
end