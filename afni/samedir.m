function issame=samedir(a,b)
% dirty check if two paths point to the same directory (as soft links may
% obscure this)
%
% ISSAME=SAMEDIR(DIR1,DIR2)
% INPUTS:
%   DIR1,DIR2    Paths to two directories
% OUTPUT:        
%   ISSAME       true iff DIR1 and DIR2 point to the same directory
%
% If either DIR1 or DIR2 is not a directory, this function returns false 
%
% Method: write a new temporary file in DIR1; if it exists in DIR2 we
% assume it's the same directory.
%
% NNO Dec 2010


if ~isdir(a), issame=false; return; end
if ~isdir(b), issame=false; return; end

if isempty(a), a='.'; end
if isempty(b), b='.'; end

a=[a filesep];
b=[b filesep];

tmpfnpat='__tmp_%d';
k=0;
while true
    tmpfn=sprintf(tmpfnpat,k);
    if ~exist([a tmpfn],'file') && ~exist([b tmpfn],'file')
        break;
    end
    k=k+1;
end

fid=fopen([a tmpfn],'w');
fprintf(fid,'foo');
fclose(fid);

issame=exist([b tmpfn],'file')>0;

delete([a tmpfn]);