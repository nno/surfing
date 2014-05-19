function theta=surfing_node_angle(v,i,j,k)

if size(v,2)~=3
    error('coordinates must be Px3');
end

coordinate_wise=nargin<=2 && size(v,2)==3 && size(i,2)==3;
index_wise=nargin>=3 && sum(size(i)~=1)<=1;

if ~(coordinate_wise || index_wise)
    error('illegal input');
end

if coordinate_wise
    theta=angle(v,i);
else
    if nargin<4
        nv=size(v,1);
        k=1:nv;
    end
    p=v(i,:);
    q=v(j,:);
    r=v(k,:);
    
    theta=angle(bsxfun(@minus,p,q),bsxfun(@minus,r,q));
end
    
    


function theta=angle(a,b)
% helper function to compute angle between two vectors
% a and b must be Px3 coordinates, theta is Px1

%theta=acos( sum(a.*b,2)./sqrt (sum(a.^2,2).*sum(b.^2,2)));
theta=acos( sum(bsxfun(@times,a,b),2)./sqrt (sum(a.^2,2).*sum(b.^2,2)));

