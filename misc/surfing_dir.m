function [d,n,p,oldd]=surfing_dir(pat)
% directory listing with support for '?' and '*' wildcard characters
%
% [D,N,OLDD]=PHOEBE_DIR(PAT)
% INPUT:
%  PAT     File name pattern: either a directory name without wildcard
%          characters, or such a directory name followed by a file name
%          thay may contain wildcard characters ('?' for any single 
%          character, and/or '*' for zero or more of any of characters).
%          PAT may also be a cell with patterns, in which case all files 
%          matching any pattern is returned, possibly multiple times if 
%          the different patterns are not mutually exclusive
% OUTPUTS:
%  D       Nx1 struct with strings of filenames with path, if N files match 
%          the pattern PAT
%  N       number of files found
%  P       path in which the files were found
%  OLDD    The output that the built-in function DIR would give you, only 
%          containing files matching the pattern PAT
%
% NNO Sep 2010
%
% See also: DIR

me=str2func(mfilename()); % immune to renaming

if iscell(pat)
    % recursion
    d=[];
    oldd=[];
    n=0;
    for k=1:numel(pat);
        [dk,nk,p,olddk]=me(pat{k});
        d=[d; dk];
        oldd=[oldd; olddk];
        n=n+nk;
    end
    return;
end

if isdir(pat)
    pat=fullfile(pat,'*');
end

% replace '?' by '*' so that builtin DIR understands 
pstar=regexprep(pat,'?','*');

% because multiple consecutive occurences of '*' makes DIR very slow, 
% replace those by a single '*'.
while strfind(pstar,'**')
    pstar=strrep(pstar,'**','*');
end

dall=dir(pstar);
n=numel(dall);

% remove leading path
[p,nm,ext]=fileparts(pat);
if numel(p)>0
    p=[p filesep()]; % add file seperator
end

% make proper regular expression. '^' and '$' mean start and end of string
preg=['^' regexptranslate('wildcard',[nm ext]) '$'];

msk=false(n,1);
for k=1:n
    fn=dall(k).name;
    msk(k)=numel(regexp(fn,preg))>0;
end

oldd=dall(msk);
n=sum(msk);
d=cell(n,1);
imsk=find(msk); % indices
for k=1:n
    d{k}=fullfile(p,dall(imsk(k)).name);
end