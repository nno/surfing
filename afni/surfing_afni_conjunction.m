function y=surfing_afni_conjunction(x,dim)
% Computes the conjunction along a dimension
%
% Y=SURFING_AFNI_CONJUNCTION(X,DIM)
% INPUTS:
%   X      P1xP2x...xPn array
%   DIM    Dimension along which conjunction is computed
% OUTPUT
%   Y      Q1xQ2x...xQn array, where Qi=1 if i==DIM and Qi=Pi otherwise
%          A nonzero value in Y means that all values along the
%          corresponding dimension in X were of the same sign. The value 
%          itself is the minimum of the absolute values along that dimension
% 
% NNO Mar 2011

sz=size(x,dim);

ax=abs(x);      % absolute value
mx=min(ax,[],dim); % minimum of absolute value
sx=(x>0)-(x<0); % sign

ssx=sum(sx,dim); % sum of signs
m=(ssx==sz)-(ssx==-sz); % -1 or 1 if all the same sign, 0 otherwise

y=bsxfun(@times,mx,m);