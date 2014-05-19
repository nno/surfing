function pos=surfing_inds2subs(siz, is)
% Multiple subscripts from linear indices; a generalization of IND2SUB
%
% P=SURFING_INDS2SUBS(SIZ,IS) returns the subscripts POS for each of the 
% linear indices in IS based on a matrix with dimensions SIZ. 
% If SIZ=[s1, ..., sn] refers to a matrix of size s1 x ... x sN, and IS is 
% a M x 1 vector that contains the linear indices, then POS is a M x N 
% matrix where each row contains the corresponding subscripts.
%
% For every index in IS that is out of bounds, the corresponding row 
% elements in P contains NaN. If IS is of type int32, however, the
% corresponding elements are set to zero.
%
% If IS is a single number (i.e. M==1) and [P1 ... PN]=IND2SUB(SIZ,IS),
% then it holds that isequal(SURFING_INDS2SUBS(SIZ,IS),[P1,...,PN])
%
% NNO June 2009, updated June 2010


rowcount=numel(is);
if size(is,2) ~= 1 
    if size(is,1) == 1
        is=is';
    else
    	error('Linear indices should be M x 1 vector');
    end
end

if ~isequal(is,floor(is)) || ~isequal(siz,floor(siz))
    error('Only integer indices and dimensions are supported');
end

if min(siz)<1
    error('Dimensions should be positive');
end

dimcount=numel(siz);

% check both inputs are same class (tested for double and int32)
clis=class(is);
clsiz=class(siz);
if ~strcmp(clis,clsiz)
    siz=cast(siz,clis);
end
dimprod=cast(prod(double(siz)),clis);

alot=1e6; % in case we have lots of indices, split them up
if rowcount>alot
    % divide in different chunks
    % rationale: REPMAT (used below) is unable to return an int32 directly, 
    % so we would run out of memory otherwise
    pos=zeros(rowcount,dimcount,clis);
    for k=1:alot:rowcount
        idxs=k:(min(k+alot-1,rowcount));
        pos(idxs,:)=surfing_inds2subs(siz,is(idxs,:));
    end
    return
end

% multiplication factors for each dimension
mply=zeros(1,dimcount,clis);
mply(1)=1;
for k=1:(dimcount-1)
    mply(k+1)=mply(k)*siz(k);
end

% see which ones are out of bounds and which ones are not
outofbounds = (is<1) | is>dimprod;
outofboundscount=sum(outofbounds);
inboundscount=rowcount-outofboundscount;

% find remainders for each of the dimensions

% NNO June 2010: this is the tricky (and hacky) part, as we want 'good' behaviour for 
% both int32 (where results are rounded) and doubles (where results are not rounded). 
rem=round((2*repmat(is(~outofbounds)-1,1,dimcount)-repmat(mply,inboundscount,1)) ./ (2*repmat(mply,inboundscount,1)));
clear diff;

% take modulo
sizrep=repmat(siz,inboundscount,1);
mods=rem - round((2*rem-sizrep+1) ./ (2*sizrep)) .* sizrep;
clear sizrep;
clear rem;
% set positions
pos=ones(rowcount,dimcount,clis);
pos(~outofbounds,:)=mods+1;
pos( outofbounds,:)=repmat(NaN,outofboundscount,dimcount);

% hack to fix 1 as input
pos(is==1,:)=repmat(ones(1,1,clis),sum(is==1),dimcount);
