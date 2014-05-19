function [p,n,v,e,e2]=afni_fileparts(fn)
% Returns the file parts of an AFNI volume file
%
% [PATH,NAME,VIEW,EXT,EXT2]=AFNI_FILEPARTS(FN) 
%
% INPUT:
%   FN      AFNI file name, optionally with path, of the form 
%           PATH '/' NAME VIEW [EXT], with VIEW in '+tlrc','+acpc','+orig'.
%           EXT can be .HEAD or .BRIK
% OUTPUTS:  
%   PATH,NAME,VIEW,EXT: As described above
%   EXT2              : 'other extension' of EXT; EXT2==.BRIK if EXT==.HEAD
%                       and vice verse
%
% NNO Dec 2010

[p,n,e]=fileparts(fn);

nm=[n e]; % file name
pluspos=find(nm=='+',1,'last');
if numel(pluspos) == 0
    if strmatch(e,{'.nii','.1D'})
        v='';
        e2='';
        return;
    else
        error('Cannot find view in %s', fn);
    end
end

n=nm(1:(pluspos-1)); % everything before '+'
ve=nm(pluspos:end); % view and extension

if numel(ve) < 5
    error('Illegal view/extension ''%s'' in %s', ve, fn);
end

idx=strmatch(ve(1:5),{'+orig','+tlrc','+acpc'});
if numel(idx)~=1
    error('View %s not understood in %s', ve(1:5),fn);
end
v=ve(1:5);

e2='';
switch(ve(6:end))
    case '.HEAD'
        e2='.BRIK';
    case '.BRIK'
        e2='.HEAD';
end

if isempty(e2)
    e='.HEAD';
    e2='.BRIK';
else
    e=ve(6:end);
end
