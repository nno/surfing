function [thetas,node_pos,v2f]=surfing_node_internal_angles(v,f)
% compute internal angles of faces for each node
%
% [thetas,node_pos,v2f]=surfing_node_angles(v,f)
%
% Inputs:
%   v              Px3 node coordinates
%   f              Qx3 face indices
%
% Outputs:
%   thetas         PxM angle (in radians), if each node is contained in at
%                  most M triangles.
%   node_pos       PxM indices in the range 1..3, indicating for each angle
%                  which column in f contains that node.
%   v2f            PxM mapping from node indices to face indices
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
%    node b made by the edges a-b and b-c, then node 1 has 6 internal
%    angles (as it has 6 faces, namely 2-1-3, 3-1-4, 4-1-5, 5-1-6, 6-1-7,
%    and 7-1-2); all other nodes have 2 internal angles.
%    * v2f(i,:) contains the indices of faces that contain node i; each
%      entry for i>1 has two values (two faces) and zeros elsewhere
%    * node_pos(i,:) contains the corresponding position in f for node i.
%      For example, node 6 is contained in faces IV and V thus v2f(6,:)
%      starts with 4(=IV) and 5(=V). f(4,:)==[1 5 6], thus the position of
%      node 6 in face 4 is the 3rd column; f(5,:)=[1 6 7], thus the
%      position of node 6 in face 5 is the 2nd column. Hence node_pos(6,:)
%      starts with [3 2] and is zero elsewhere.
%    * theta contains the corresponding angles. For example the angles for
%      node 6, 1-5-6 and 1-6-7, are both around 1.11 radians, or
%      about 63 degrees.
%
%    Using this function gives
%    >> [thetas,node_pos,v2f]=surfing_node_internal_angles(v,f)
%
%    thetas =
%
%     1.1071    1.1071    0.9273    1.1071    1.1071    0.9273
%     0.9273    1.1071         0         0         0         0
%     1.1071    1.1071         0         0         0         0
%     0.9273    1.1071         0         0         0         0
%     1.1071    0.9273         0         0         0         0
%     1.1071    1.1071         0         0         0         0
%     0.9273    1.1071         0         0         0         0
%
%
%    node_pos =
%
%      1     1     1     1     1     1
%      2     3     0     0     0     0
%      3     2     0     0     0     0
%      3     2     0     0     0     0
%      3     2     0     0     0     0
%      3     2     0     0     0     0
%      3     2     0     0     0     0
%
%
%    v2f =
%
%      1     2     3     4     5     6
%      1     6     0     0     0     0
%      1     2     0     0     0     0
%      2     3     0     0     0     0
%      3     4     0     0     0     0
%      4     5     0     0     0     0
%      5     6     0     0     0     0
%
% Notes:
%   - If node K is contained in J faces, then:
%     * each of the three outputs has non-zero entries at (K,1:J) and zero
%       entries at (K,(J+1,end))
%     * v2f(K,I)==G means that face G contains node K, and more
%       specifically, f(G,node_pos(K,I))==I, with the angle at node I in
%       face G being angles(K,I)
%   - To convert angles to degrees, multiply by (180/pi).
%
% NNO Apr 2014

% mapping: node to face
v2f=surfing_nodeidxs2faceidxs(f');

% mapping: face to angle
f2a=surfing_face2angle(v,f);

nv=size(v,1);
node_idxs=(1:nv)';
nnbrs_max=size(v2f,2);

% temporary variables
angles=zeros(nv,3);
fs=zeros(nv,3);

% space for output
thetas=zeros(nv,nnbrs_max);
node_pos=zeros(nv,nnbrs_max);
for col=1:nnbrs_max
    % for each node, get a face containing that node
    f_idxs=v2f(:,col); % indices of faces that contain the node

    % apply mask for nodes with less than 'col' neighbors
    msk=f_idxs>0;

    % angles associated with each node
    angles(:)=0;
    angles(msk,:)=f2a(f_idxs(msk),:);

    fs(:)=0;
    fs(msk,:)=f(f_idxs(msk),:); % nodes in the face

    for j=1:3
        m=fs(:,j)==node_idxs & msk;
        assert(all(thetas(m,col)==0));
        thetas(m,col)=angles(m,j);
        node_pos(m,col)=j;
    end
end
