function err=surfing_check_surface(v,f,raise_if_error,level)
% perform basic checks for surface
%
% err=surfing_check_surface(v,f,raise_if_error)
%
% Inputs:
%  v                Px3 node coordinates
%  f                Qx3 faces containing node indices
%  raise_if_error   If true (the default) then an error is raised if
%                   the surface is not proper. This can occur when
%                   - v or f are of the wrong shape
%                   - indices in f are inconsistent with v, or non-finite
%  level            If >=2, then it is also verified that duplicate edges
%                   do not exist.
% Returns
%  err              The error message if a surface is not proper, or empty
%                   otherwise
%
% NNO Apr 2014


if nargin<4 || isempty(level), level=1; end
if nargin<3 || isempty(raise_if_error), raise_if_error=true; end


err='';

while true
    [nv,three]=size(v);
    if three~=3
        err='vertices must be Px3';
        break;
    end

    [nf,three]=size(f);
    if three~=3
        err='faces must be Qx3';
        break;
    end

    f_arr=f(:);
    if min(f_arr)<1
        err=sprintf('min face index %d < 1', min(f_arr));
        break;
    end

    if max(f_arr)>nv
        err=sprintf('max face index %d > %d', max(f_arr), nv);
        break;
    end

    if ~all(isfinite(f_arr))
        err=sprintf('face indices must be finite');
        break;
    end

    if ~isequal(round(f_arr),f_arr)
        err=sprintf('face indices must be integers');
        break;
    end

    v_none=true(nv,1);
    v_none(f(:))=false;
    if any(v_none)
        idx=find(v_none,1);
        err=sprintf('%d nodes (including node %d) not indexed by faces', ...
                                sum(v_none), idx);
        break;
    end

    [sf,i]=sortrows(sort(f')');
    m=all(diff(sf)==0,2);
    if any(m)
        idx=i(find(m,1));
        idx2=i(find(m,1)+1);
        err=sprintf('duplicate face %d and %d with nodes %d, %d, %d',...
                    idx,idx2,f(idx,:));
        break;
    end


    if level>1
        fs=1:nf;
        ii=zeros(nf*3,1);
        jj=zeros(nf*3,1);
        ss=zeros(nf*3,1);

        for k=1:3
            idxs=(1:nf)+(k-1)*nf;
            ii(idxs)=f(:,k);
            jj(idxs)=f(:,mod(k+1,3)+1);
            ss(idxs)=1;
        end

        m=ii>0;
        edge2count=sparse(ii(m),jj(m),ss(m),nv,nv,sum(m));
        [mxi,i]=max(edge2count);
        [mxj,j]=max(mxi);
        if mxj>1
            for k=1:3
                msk=f(:,k)==(i(j)+0) & f(:,mod(k+1,3)+1)==(j+0);
                if any(msk)
                    err=sprintf('duplicate edge (%d,%d) at faces %d and %d',...
                                j, i(j), find(msk,2));

                    break
                end
            end
            assert(~isempty(err), 'this should not happen');
        end
    end

    break;
end

if raise_if_error && ~isempty(err)
    error(err);
end

