function R=surfing_suma_surfacefiles(specfn,trgdir,cpmethod)
% Finds AFNI SUMA surface files, and optionally copies them to another
% directory. 
%
% [R]=SUMA_SURFACEFILES(SPECFN[,TRGDIR,CPMETHOD])
% INPUTS:
%   SPECFN:     SUMA .spec file; or if a directory is given, then all .spec 
%               files in that directory are processed.
%   TRGDIR:     Directory where files are copied to; use [] to not copy the
%               files
%   CPMETHOD:   'ln' or 'cp' for hard link or 'normal' copy
% OPTIONAL OUTPUT:
%   R           Struct with fields .surfs, .specfile, .shell, .anat1,
%               .anat2, .allfiles, .dir, and .(S) where S is any name of a
%               surface
%               
%               
% If TRGDIR is specified, then function copies all ASCII Freesurfer files 
% defined in the .spec file. If there is a .sh file that calls suma, then 
% the corresponding SurfVol anatomical files are copied too.
%
% Note that this function assumes certain naming conventions for the 
% different files (as specified in the code). It should work for files
% generated with SURFING_SUMA_ALIGN and SURFING_SUMA_MAKESPEC.
%
% See also: SURFING_SUMA_ALIGN, SURFING_SUMA_MAKESPEC
%
% NNO Dec 2010

me=str2func(mfilename());

if nargin<3
    cpmethod='ln';
end


if nargin<2 || isempty(trgdir) || ~ischar(trgdir)
    trgdir=NaN;
end

if ~strmatch(cpmethod,{'ln','cp'})
    error('Not supported cpmethod: %s', cpmethod);
end

[p,f,e]=fileparts(specfn);

if isdir(specfn)
    d=surfing_dir(fullfile(specfn,'*.spec'));
    n=numel(d);
    if n==0
        error('No spec files found in %s', specfn);
    end
    R=cell(n,1);
    for k=1:n
        R{k}=me(d{k}, trgdir, cpmethod);
    end
    return
end

if ~exist(specfn,'file')
    error('Spec file %s not found', specfn);
end



if ischar(trgdir) && samedir(p,trgdir)
    % note: we cannot check is specfile exists, because it could be an
    % older one that we want to overwrite. Therefore we really have to
    % check that the sourcefolder is the same as the destination folder.
    warning('Sourcedir and target are the same, will not overwrite specfile')
    cpmethod='';
elseif ischar(trgdir) && ~isdir(trgdir)
    mkdir(trgdir);
end

surffiles=cell(0);

keys={'FreeSurferSurface','LabelDset','SurfaceName','LocalDomainParent'};
ignorekeys={'SAME'};

R=struct();

pat='\s*(?<key>\S+)\s*=\s*(?<val>[/\.\S]+)\s*';
fid=fopen(specfn);
allhemis=cell(0);
while true
    line=fgetl(fid);
    if isnumeric(line)
        break;
    end

    m=regexp(line,pat,'names');
    if ~isempty(m) 
        for j=1:numel(keys)
            if isnan(m.key)
                break;
            end
            if strcmpi(keys{j},m.key) 
                for k=1:numel(ignorekeys)
                    if findstr(ignorekeys{k},m.val)
                        m.key=NaN;
                    end
                end
                if isnan(m.key) % matlab does not jumping out of multiple loops
                    break;
                end
                    
                
                surffiles{end+1}=m.val;
                mhemi=regexp(m.val,'.*([lr]h).*\.asc','tokens');
                if numel(mhemi)==1
                    hemi=mhemi{1}{1};
                    allhemis{end+1}=hemi(1);
                end
                
                mico=regexp(m.val,'\D*ico\D*(\d+)\D*','tokens');
                if numel(mico)==1
                    ico=mico{1}{1};
                    R.ico=ico;
                end
                m.key=NaN; % break out
            end
        end 
    end

end

allhemis=unique(allhemis);
R.hemi=[allhemis{:}];

surffiles=unique(surffiles);

% add to get easy access to file names
surfnames={'semiinflated','tqinflated','smoothwm','pial','inflated','intermediate'};
for k=1:numel(surffiles)
    m=regexp(surffiles{k},surfnames);
    for j=1:numel(m)
        if ~isempty(m{j})
            R.(surfnames{j})=[p '/' surffiles{k}];
            break;
        end
    end
            
end

R.surfs=surffiles;

surffiles{end+1}=[f e]; % spec file itself
R.specfile=surffiles{end};


if exist('hemi','var') && exist('ico','var')
    n=0;
    pats={[p '/*' ico '_*' hemi '*.sh'],[p '/*' hemi '*' ico '_*.sh']};
    for j=1:numel(pats)
        d=dir(pats{j});
        n=n+numel(d);
        
        if numel(d)==1
            shellfn=[p '/' d(1).name];
            break;
        end
    end
    
    if n==1
        fprintf('Found unique shell script: %s\n', shellfn);
        surffiles{end+1}=d(1).name;
        R.shell=surffiles{end};
        pat='.*suma.*\s+-sv\s+(\S+).*';
        
        fid=fopen(shellfn);
        while true
            line=fgetl(fid);
            if isnumeric(line)
                break;
            end
            
            manat=regexp(line,pat,'tokens');
            if ~isempty(manat)
                manat=manat{1}{1};
                [dummy f v e e2]=afni_fileparts(manat);
                surffiles{end+1}=[f v e];
                surffiles{end+1}=[f v e2];
                R.anat1=surffiles{end-1};
                R.anat2=surffiles{end};
                break;
            end
        end
        fclose(fid);
        
        if ~isfield(R,'anat1')
            error('Not found anat');
        end
    end
end

R.allfiles=surffiles;
R.dir=[p '/'];

if ~ischar(trgdir)
    return;
end
    
surffiles=unique(surffiles);
n=numel(surffiles);
cmd='';
if isempty(trgdir)
    trgdir='.';
end
for k=1:n
    fn=surffiles{k};
    src=[p '/' fn];
    trg=[trgdir '/' fn];
    switch cpmethod
        case 'ln'
            cmd=sprintf('%sln -f %s %s || cp -f %s %s;',cmd,src,trg,src,trg);
        case 'cp'
            cmd=sprintf('%scp -f %s %s;',src,trg);
            
    end
end

unix(cmd);


    

