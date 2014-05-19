function theta=surfing_face2angle(v,f)
% compute angles in faces of a surface
%
% theta=surfing_face2angle(v,f)
%
% Inputs:
%   v              Px3 node coordinates
%   f              Qx3 face indices
%
% Output:
%   theta          Qx3 angles (in radians). theta(i,j)=t means that node
%                  p=f(i,j) in face i makes an angle t with the other two
%                  nodes in face i.
%
% Notes:
%   - to convert the angles from radians to degrees, multiply theta by 
%     (180/pi)
%
% NNO Apr 2014

surfing_check_surface(v,f);
nv=size(v,1);
nf=size(f,1);

theta=zeros(nf,1);
for i=1:3
    j=mod(i,3)+1;
    k=mod(i+1,3)+1;
    
    p=v(f(:,i),:)-v(f(:,j),:);
    q=v(f(:,i),:)-v(f(:,k),:);

    theta(:,i)=angle(p,q);
end

function theta=angle(a,b)
% helper function to compute angle between two vectors
% a and b must be Px3 coordinates, theta is Px1

theta=acos( sum(a.*b,2)./sqrt (sum(a.^2,2).*sum(b.^2,2)));

    

