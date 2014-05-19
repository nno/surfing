function [sidxs,r2]=surfing_selectkfirstidxs(k,idxs,inclmsk)
% selects approximately the first k non-zero unique indices
%
% SIDXS=SURFING_SELECTKFIRSTIDXS(K,IDXS,[ROWHINT])
%
% INPUT:
%  K:       number of indices to select
%  IDXS:    NxP matrix with indices, corresponding to N nodes with at most P
%           indices each.
% OUTPUT:
%  SIDXS:   K2x1 matrix with unique indices, taken from IDXS. SIDXS contains
%           the indices from IDXS, starting with the first row, the second
%           row, etc, until K2 (approximately K) indices are found that 
%           are all unique and non-zero. Either all values from a row are
%           selected or none. The function tries to find the best number of
%           rows in an unbiased way. 
% ROWCOUNT: number of rows selected from IDXS.
%
% NNO May 2010

[n,p]=size(idxs);


r1=0; % lower bound for number of rows
r3=n; % upper "                      "

rowmsk=false(n,p);

totalunqcount=numel(unique([0; idxs(:)]));
if totalunqcount<k
    error('Cannot select %d indices, because only %d unique indices', k, totalunqcount-1);
end

while true
    r2=ceil((r1+r3)/2); % select up to and including r2 rows

    rowmsk(1:r2,:)=idxs(1:r2,:)>0;
    rowmsk((r2+1):end,:)=false;
    sidxs=unique(idxs(rowmsk));
    nidxs=numel(sidxs);

    if nidxs==k % exact match, we're done
        break;
    end
    
    
    % see which indices we select with one row fewer
    rowmsk(r2,:)=false;
    sidxs2=unique(idxs(rowmsk));
    nidxs2=numel(sidxs2);
    
    if nidxs2<k && k<nidxs % k falls in between, we're almost done
    
        % see if it's better to drop the last (i.e. r2-th) row
        curvs2=(k-nidxs2)-(nidxs-k); % which one

        if sum([(k-nidxs2),(nidxs-k)]>0) ~= 2
            error('This should not happen')
        end
        
        if curvs2==0 && mod(r2,2)==0 
            % it's a tie. select pseudorandomly based on whether number of
            % rows is even or not (this should be unbiased)

            curvs2=-1; % set negative to drop last (i.e. r2-th) row
        end

        if curvs2<0 && r2>1
            % drop the last row (because nidxs2 is closer to k than
            % nidxs, or because it's a tie and nidxs2 was lucky to win)
            r2=r2-1;
            sidxs=sidxs2;
        end
        
        break; % and we're done
    end

    if nidxs<k
        r1=r2;
    else
        r3=r2;
    end
end