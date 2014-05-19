function b=surfing_mask_borders(msk,nbrs)
% returns nodes bordering a masked region
%
% NNO Aug 2011


[nv,nn]=size(nbrs);
if iscell(msk) 
    msk=unique([msk{:}]);
end

if isnumeric(msk) 
    if isequal(round(msk),msk) && min(msk)>=1 && max(msk) <= nv
        mskidxs=msk;
        msk=false(nv,1);
        msk(mskidxs)=true;
    else
        msk=logical(msk);
    end
elseif ~islogical(msk)
    error('Illegal mask type');
end

if numel(msk) ~= nv, error('Size mismatch between mask and neighbours'); end
    
nbrs=[ones(1,nn); 1+nbrs]; % add dummy 
nbrsh=sum(nbrs>1,2); % number of neighbours for each value

mskn=[false; msk(:)]+0; % convert to numeric; incorporate dummy
mskh=zeros(nv+1,1);

for k=1:nn
    mskh=mskh+mskn(nbrs(:,k));
end

b=mskh(2:end)<nbrsh(2:end) & msk;




        
    
    
    