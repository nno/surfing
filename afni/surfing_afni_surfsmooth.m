function [cmd,fnout]=surfing_afni_surfsmooth(surfdatapat,specfn,targetfwhm,varargin)
% Wrapper function to run AFNI SUMA SurfSmooth
%
% CMD=SURFING_AFNI_SURFSMOOTH(SPAT,SPECFN,FWHM,...)
% INPUTS:
%   VPAT      Pattern of .niml.dset surface files to smooth (wildcards allowed)
%   SPECFN    SUMA spec file, or a directory with such a file. In the
%             latter case the function makes an educated guess about
%             hemispheres
%   FWHM      Target smoothness expressed in FWHM. If FWHM is negative,
%             then the data is smoothed with -FHWM irrespective of the
%             original smoothness
%   ...       Any option for SurfSmooth, either as struct or key-value,
%             using syntax supported by SURFING_STRUCT.
%             Defaults are smoothing method 'HEAT_07', master is detrended,
%             and the intermediate surface is used for smoothing
%
% See also: SURFING_STRUCT, SURFING_AFNI_RUNBINARY,
% SURFING_SUMA_SURFACEFILES
%
% NNO Mar 2011  


% set defaults
Df=struct();
Df.met='HEAT_07';
Df.detrend_master=true;
Df.surf_A='intermediate';
Df.replexp_={'.niml.dset','_sblur%d.niml.dset'};
me=str2func(mfilename()); % make immune to renaming

if iscell(specfn)
    cmd=cell(size(specfn));
    for k=1:numel(specfn)
        cmd{k}=me(surfdatapat,specfn{k},targetfwhm,varargin{:});
    end
    return;
elseif isstruct(specfn)
    specfn=[specfn.dir specfn.specfile];
elseif ischar(specfn) && isdir(specfn)
    specfn=surfing_suma_surfacefiles(specfn);
    cmd=me(surfdatapat,specfn,targetfwhm,varargin{:});
    return;
end

R=surfing_suma_surfacefiles(specfn);
[fns,n,surfdir]=surfing_dir(surfdatapat);

cmd=sprintf('cd %s', surfdir);

if n==0
    warning('No files found matching %s\n', surfdatapat);
    return
end

if targetfwhm<0
    Df.fwhm=-targetfwhm;
else
    Df.target_fwhm=targetfwhm;
end
if round(targetfwhm)==targetfwhm
    blurstr=sprintf('%d',round(targetfwhm));
else
    blurstr=sprintf('%.1f',targetfwhm);
end

fnout=cell(0);

for j=1:n
    fullfn=fns{j};
    [p,nm,ext]=fileparts(fullfn);
    
    if isempty(strfind(nm,[R.hemi 'h']))
        msg=sprintf('Skipping %s: did not find expected hemisphere pattern %sh', fullfn, R.hemi);
        warning(msg);
        cmd=sprintf('%s; echo %s',cmd,msg);
        continue
    end
    
    % fileparts does not work with double extensions
    trgext='.niml.dset';
    trgextstart=strfind([nm ext],trgext);
    if isempty(trgextstart)
        error('Expected file with extension %s in %s, but not found',trgext,fullfn);
    end
 
    Df.input=[nm ext];
    C=surfing_struct(Df,varargin{:});
    
    replwith=regexprep(C.replexp_{2},'%d',blurstr);
    outputfn=regexprep([nm ext],C.replexp_{1},replwith);

    Df.output=outputfn;
    if targetfwhm>0
        Df.blurmaster=C.input;
    end
    Df.spec=R.specfile;
    Df.overwrite=true;
    
    C=surfing_struct(Df,varargin); % join all options
    opt=surfing_afni_opts2string(C); % convert to string
    
    cmd=sprintf('%s; echo; echo Blurring %d / %d for %sh: %s using %s',cmd,j,n,R.hemi,C.input,C.spec);
    cmd=sprintf('%s;%s %s',cmd,surfing_afni_runbinary('SurfSmooth'),opt);
    cmd=sprintf('%s; echo Completed %d / %d: %s; echo',cmd,j,n,C.input);
    
    fnout{end+1}=fullfile(p,outputfn);
end

unix(cmd,'-echo');