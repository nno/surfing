function is=surfing_subs2inds(siz, pos)
% Linear indices from multiple subscripts; a generalization of sub2ind
%
% IS=SURFING_SUBS2INDS(SIZ,P) returns the linear indices for each row of 
% subscripts in P based on a matrix with dimensions SIZ. 
% If SIZ=[s1, ..., sn] refers to a matrix of size s1 x ... x sN, and P is 
% an M x N where each row contains the subindices referring to a single 
% element, then IS is a M x 1 vector with linear indices.
%
% For indices in P that are out of bounds, the corresponding element in IS
% is NaN. If P is of type int32, however, the corresponding element is set 
% to zero.
%
% If POS=[p1 ... pN] (i.e. M==1) and [p1,...,pN]=sub2ind(SIZ,IS), then
%    [p1,...,pN]=SURFING_SUBS2INDS(SIZ,IS)
%
% NNO May 2009, updated June 2010

dimcount=numel(siz);
[poscount,dimcount2]=size(pos);

if dimcount ~= dimcount2
    error('Number of dimensions do not match');
end

% check both inputs are same class (tested for double and int32)
clpos=class(pos);
clsiz=class(siz);
if ~strcmp(clpos,clsiz)
    siz=cast(siz,clpos);
end

% multiplication factors for the different positions
mply=zeros(1,dimcount,clpos);
mply(1)=1;
for k=1:(dimcount-1)
    mply(k+1)=mply(k)*siz(k);
end

alot=1e6;
if poscount>alot
    is=zeros(poscount,1,clpos);
    for k=1:alot:poscount
        idxs=k:(min(k+alot-1,poscount));
        is(idxs)=surfing_subs2inds(siz,pos(idxs,:));
    end
    return
end

is=sum((pos-1) .* repmat(mply,poscount,1),2,'native')+1;

beyond=repmat(siz,poscount,1)-pos;
outofbounds=sum((pos<1)|(beyond<0),2)>0;
is(outofbounds)=NaN;