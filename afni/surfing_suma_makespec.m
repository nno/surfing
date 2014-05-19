function cmd=surfing_suma_makespec(varargin)
% Makes a AFNI SUMA spec file and script to start SUMA
%
% SURFING_SUMA_MAKESPEC('aligndir',D,surfvolalfn',S,'icold',I)
%
% INPUTS:
%   D      Directory that contains surface files
%   S      Anatomical file (e.g. 's01_SurfVol_fix+orig')
%   I      MapIcosehedron number of subdivisions
%
% This script supports surfaces for 'smoothwm','intermediate','pial',
% 'semiinflated','tqinflated','inflated'.
%
% As of Mar 2011, it also accepts the output from SURFING_SUMA_ALIGN
%
% NNO Jan 2010

D=struct();
D.icopat='ico%d_';
C=surfing_struct(D,varargin,'?',{'aligndir','surfvolalfn','icold','icopat'});
surfing_struct(C);

[p,n,view]=afni_fileparts(surfvolalfn);

if ~isempty(aligndir)
    aligndir=[aligndir '/'];
end

switch view
    case {'orig','+orig'}
        pat='_al.asc';
    case {'tlrc','+tlrc'}
        pat='_tlrc.asc';
    otherwise
        error('view should be orig or tlrc, but found %s', view);
end

view=regexprep(view,'+',''); % remove '+' sign, if present

ifxs={'smoothwm','intermediate','pial','semiinflated','tqinflated','inflated'}; % surface infixes

cmd=surfing_afni_runbinary();
for hemi='lr'
    fnroot=sprintf(['%sh_' icopat '%s'],hemi,icold,view);
    headerlines=sprintf('# Created %s\nGroup = all\n',datestr(now));
    lines=[];
    
    n=numel(ifxs);
    for k=1:n
        ifx=ifxs{k};
        surffn=sprintf([icopat  '%sh.%s%s'],icold,hemi,ifx,pat);
        if ~exist([aligndir surffn],'file')
            surfpat=sprintf([aligndir icopat '%sh.%s*%s'],icold,hemi,ifx,pat);
            surfdir=surfing_dir(surfpat);
            if numel(surfdir)==1
                surffn=surfdir{1};
            else
                warning('pattern %s not in %s',surfpat,aligndir)
                continue;
            end
        end
        
        if k==1
            ldp='SAME';
            parent=surffn;
        else
            ldp=parent;
        end
        
        
        lines=sprintf(['%s\nNewSurface\nSurfaceFormat = ASCII\nSurfaceType = FreeSurfer\n' ...
					'FreeSurferSurface = %s\nLocalDomainParent = %s\n'...
					'SurfaceState = %s\nEmbedDimension = 3\n\n'],...
                    lines,surffn,ldp,ifx);
        
        headerlines=sprintf('%sStateDef = %s\n',headerlines,ifx);
    end
    
    lines=sprintf('%s\n%s',headerlines,lines);

    
    specfn=sprintf('%sspec_%s.spec',aligndir,fnroot);
    seesumafn=sprintf('%s%s_seesuma.sh',aligndir,fnroot);
    
    seesumalines=sprintf('export SUMA_AllowDsetReplacement=YES;killall afni\nafni -niml &\nsuma -spec %s -sv %s', nopath(specfn), nopath(surfvolalfn));
    
    writestr(specfn,lines);
    writestr(seesumafn,seesumalines);
    
    cmd=sprintf('%schmod u+x %s;',cmd,seesumafn);
    unix(cmd);
    
    fprintf(['\nSpec file for %sh written to %s\nTo view in SUMA and AFNI, start a shell and run:\n' ...
        '  cd %s\n  ./%s\n\n'], hemi, specfn, aligndir, nopath(seesumafn));
end

function ne=nopath(fn)
[p,n,e]=fileparts(fn);
ne=[n e];

function writestr(fn,s)
fid=fopen(fn,'w');
fwrite(fid,s);
fclose(fid);
