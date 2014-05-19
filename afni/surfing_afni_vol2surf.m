function cmd=surfing_afni_vol2surf(volpat,specfn,varargin)
% Wrapper function to run AFNI SUMA Vol2Surf
%
% CMD=SURFING_AFNI_VOL2SURF(VPAT,SPECFN,...)
% INPUTS:
%   VPAT      Pattern of HEAD volumefiles to project (wildcards allowed)
%   SPECFN    SUMA spec file, or a directory with such a file. In the
%             latter case the function makes an educated guess about
%             hemispheres
%   ...       Any option for 3dVol2Surf, either as struct or key-value,
%             using syntax supported by SURFING_STRUCT.
%             Defaults are the 'ave' map function, with ten steps from pial
%             to smoothwm surface
%
% See also: SURFING_STRUCT, SURFING_AFNI_RUNBINARY,
% SURFING_SUMA_SURFACEFILES
%
% NNO Mar 2011

% set defaults
Df=struct();
Df.mapfunc='ave';
Df.f_steps=10;
Df.f_index='nodes';
Df.surf_A='pial';
Df.surf_B='smoothwm';

me=str2func(mfilename()); % make immune to renaming

if iscell(specfn)
    cmd=cell(size(specfn));
    for k=1:numel(specfn)
        cmd{k}=me(volpat,specfn{k},varargin{:});
    end
    return;
elseif isstruct(specfn)
    specfn=[specfn.dir specfn.specfile];
elseif ischar(specfn) && isdir(specfn)
    specfn=surfing_suma_surfacefiles(specfn);
    cmd=me(volpat,specfn,varargin{:});
    return;
end

R=surfing_suma_surfacefiles(specfn); % get surface files

[fns,n,voldir]=surfing_dir(volpat);

%cmd=sprintf('cd %s', voldir);
cmd=sprintf('cd %s', R.dir);

if n==0
    warning('No files found matching %s\n', volpat);
    return
end

for j=1:n
    fullfn=fns{j};
    [p,nm,vw]=afni_fileparts(fullfn);
    
    C=struct();
    C.out_niml=[p '/' R.hemi 'h_' nm '.niml.dset'];
    C.sv=R.anat1;
    C.grid_parent=[p '/' nm vw];
    C.spec=R.specfile;
    
    C=surfing_struct(Df,C,varargin); % join all options
    C=makelast(C,{'spec','sv','grid_parent'}); % some fields as last args
    opt=surfing_afni_opts2string(C); % convert to string
    
    cmd=sprintf('%s;rm %s; echo; echo "Mapping %d / %d for %sh using %s: %s -> %s"',cmd,C.out_niml,j,n,R.hemi,C.spec,C.grid_parent,C.out_niml);
    cmd=sprintf('%s;%s %s',cmd,surfing_afni_runbinary('3dVol2Surf'),opt);
    cmd=sprintf('%s; echo "Completed %d / %d: %s"; echo',cmd,j,n,C.sv);
end

unix(cmd,'-echo');
surfing_suma_surfacefiles(specfn,voldir);

function C=makelast(C,fs) % fs are fieldnames that are set last in C

n=numel(fs);
for k=1:n
    f=fs{k};
    v=C.(f);
    C=rmfield(C,f);
    C.(f)=v;
end

