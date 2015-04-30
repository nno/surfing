function [nbrs,dst]=surfing_surface_nbrs(f,v,sort_)
% Finds the neighbours on a surface
%
% [nbrs[,dst]]=surfing_surface_nbrs(f[,v[,sort_]])
%
% INPUT:
%  F      Nx3 face matrix for N faces
%  V      optional Mx3 coordinate matrix for M vertices. Only required if
%         distances are to be returned
%  SORT_  If false (the default) the distances returned are in arbitrary
%         order. If true then the distances are sorted in each row.
%
% OUTPUTS:
%  NBRS   NxP neighbour matrix, if each node has at most D neighbours
%         If node K has Q neighbors, then the K-th row contains the
%         Q indices with the neighbors of node K followed by (P-Q) zeros.
%  DST    NxP distance matrix containing the euclidian distances to each
%         neighboring node, padded with zeros as in NBRS.
%
%
% NNO Sep 2011

[n,k]=size(f);
if n==3 && k~=3
    warning('surfing:surface_nbrs','surfing_nbrs expected Nx3 face matrix; will transpose ');
    f=f';
    [n,k]=size(f);
end

[fi,nfi]=surfing_invertmapping(f); % vertices to facs
[nv,mx]=size(fi); % number of vertices, and max number of neighbours per face
allnbrs=zeros(nv,mx*k); % vertices contained in the faces containing the center vertex
for j=1:mx
    msk=nfi>=j; % omit zero values
    allnbrs(msk,(j-1)*k+(1:k))=f(fi(msk,j),:);
end

unqnbrs=surfing_uniqueidxsperrow(allnbrs);
nbrs=zeros(size(unqnbrs)-[0 1]);

return_distances=nargout>=2;
if return_distances
    if nargin<2
        error('need vertices to compute distances');
    end
    if nargin<3
        sort_=false;
    end
    nv=size(v,1);
    if nv ~= size(nbrs,1)
        % check size
        error('nbrs (%d) and vertex count (%d) mismatch',size(nbrs,1),nv);
    end
    % allocate space for distances
    dst=zeros(size(nbrs));

    % skip nodes with non-finite values
    skip_node_mask=any(~isfinite(v),2);

    % set coordinates to NaN
    v(skip_node_mask,:)=NaN;
end

for j=1:nv
    % remove self and empty values
    msk=unqnbrs(j,:) ~= 0 & unqnbrs(j,:) ~= j;
    nbrs(j,1:sum(msk))=unqnbrs(j,msk);
end

if return_distances
    n_max=size(nbrs,2);

    for k=1:n_max
        % consider nodes with exactly k neighbors
        nbr=nbrs(:,k);
        m=nbr>0;

        % compute euclidian distance
        delta=v(m,:)-v(nbr(m),:);
        dst(m,k)=sum(delta.^2,2).^.5;
    end

    % remove nodes with NaN
    remove_nodes_mask=isnan(dst);
    nbrs(remove_nodes_mask)=0;


    if sort_
        % count number of neighbors for each node
        count=sum(nbrs>0,2);

        % space for sorted distances
        sorted_dst=zeros(size(dst));
        sorted_nbrs=zeros(size(nbrs));

        for k=1:n_max
            % consider nodes with exactly k neighbors
            m=find(count==k);
            if isempty(m)
                continue;
            end

            % get distances and neighbors
            d=dst(m,1:k);
            nbr=nbrs(m,1:k);

            % sort distances row-wise
            [unused,i]=sort(d',1);

            % because the sorted indices are set per row they
            % cannot be indexed directly. Rather have a double for loop
            % to find index p in column j
            for j=1:k
                for p=1:k
                    mj=i(j,:)==p;
                    mmj=m(mj);
                    sorted_dst(mmj,j)=d(mj,p);
                    sorted_nbrs(mmj,j)=nbr(mj,p);
                end
            end
        end

        nbrs=sorted_nbrs;
        dst=sorted_dst;
    end


end







