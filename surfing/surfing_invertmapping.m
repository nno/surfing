function [y2x,n]=surfing_invertmapping2(x2y,ymask)
% inverts a (possibly non-injective) mapping from N1 to N1, where N1 is the
% set of positive integers.
%
% [Y2X,NX]=SURFING_INVERTMAPPING(X2Y) 
% INPUT:
%   X2Y:   NxK vector with postive (>=1) integer values
%   YMASK: Optional S-element vector with mask of values in X2Y to include
%          K=X2Y(I,J) is excluded if K<=S and YMASK(K) is logically false.
%          
% OUTPUTS:
%   Y2X:   PxQ matrix with postive integer values, where P==MAX(X2Y(:)), and 
%          Q is the maximum number of occurences of any value in X2Y.
%   NX:    Px1 matrix with the number of occurences of elements in X2Y.
%
% Properties of X2X: if Y2X(I,J)==K, then either:
%  - K<=0, which means that X2Y contains less than J occurences of I.
%  - K>0 , which means that X2Y(K,M)==I for certain M.
%  
% Example: SURFING_INVERTMAPPING([1 6 4 1 3 4]') returns [1 4
%                                                         0 0
%                                                         5 0
%                                                         3 6
%                                                         0 0
%                                                         2 0]
%
% NNO Nov 2009, May 2010
% updated Sep 2010 to support multiple columns.

if ~isequal(floor(x2y),x2y) || sum(isinf(x2y(:)))>0
    error('Only integers supported');
end

[p,q]=size(x2y);
x2ymax=max(x2y(:)); % maximum value in x2y
n=histc(x2y(:),1:x2ymax); % number of occurences of each value
maxn=max(n); % max number of occurences, i.e. number of columns in output

% set y mask, if not given
if nargin<2 || isempty(ymask)
    ymask=true(x2ymax,1);
elseif  ~isequal(floor(ymask),ymask) || sum(isinf(ymask(:)))>0
    error('Only integers for ymask supported');
end
nymask=numel(ymask);
    
% allocate space for output
y2x=zeros(x2ymax,maxn);


% keep track of first free position in each row
freepos=zeros(x2ymax,1);

for k=1:p 
    for j=1:q
        xval=x2y(k,j);
        if xval > 0 && (xval > nymask || ymask(xval))
            freepos(xval)=freepos(xval)+1;
            y2x(xval,freepos(xval))=k;
        end
    end
end
