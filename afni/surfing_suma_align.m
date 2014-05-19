function [cmd,C]=surfing_suma_align(varargin)
% Alignment of Freesurfer surfaces to an anatomical volume using AFNI
% 
% Necessary parameters:
%   action:     'tosuma', 'mapico', 'moresurfs', 'alignoriganat',
%              'applyoriganat', 'all', or ''
%   origvolfn:  path and filename of original anatomical volume that must be
%              aligned to the functional data.
%   fssurfdir:  Freesurfer "surf" dir 
%   subjid:     subject id
%   aligndir:   directory to do alignment in (and put results in)
% These parameters should be given as in SURFING_STRUCT.
%
% This function will execute several AFNI and SUMA binaries to convert the
% Freesurfer surfaces to AFNI/SUMA format (the 'tosuma' step), convert
% these to a common topology ('mapico'), generate additional surfaces such
% as intermediate and semi-inflated ('moresurfs'), find the alignment
% between the Freesurfer surfaces to the functional data ('alignoriganat'),
% and apply this alignment to the surfaces and check the algignment
% ('applyoriganat').
%
% Note: this function is *experimental* and assumes that files are stored
% with certain names specified in the body of this function.
%
% See also SURFING_STRUCT, SURFING_SUMA_MAKESPEC.
%
% NNO Mar 2011 

me=str2func(mfilename());

C=surfing_struct(varargin{:},'?',{'fssurfdir','origvolfn','subjid','aligndir','action'});

matI='"MATRIX(1,0,0,0,0,1,0,0,0,0,1,0)"'; % identity matrix
matlpirai='"MATRIX(-1,0,0,0,0,-1,0,0,0,0,1,0)"'; % swap between lpi and rai

D=struct();
D.action='action not defined';
D.afnidir='/sw/afni';
D.shellinit='. ~/.bash_profile;export OMP_NUM_THREADS=6;';
D.hemis='lr';
D.moresurfnames={'smoothwm','pial','intermediate';'intermediate','inflated','semiinflated';'semiinflated','inflated','tqinflated'};
D.allsurfnames=unique(D.moresurfnames); % assumes all names have been used in moresurfnames

C=surfing_struct(D,varargin{:});
% set some variables that are used below
D.sumasubjid=C.subjid;
D.icopat='ico%d_';
C=surfing_struct(D,varargin{:});

D.icold=100;
D.hemis='lr';
D.sumadir=[C.fssurfdir '/SUMA/'];
D.matfile=sprintf('%s_SurfVol_al_mat.aff12.1D',C.sumasubjid);
D.matfileinv=sprintf('%s_SurfVol_al_mat_inv.aff12.1D',C.sumasubjid);
D.matfile3x4=['3x4_' D.matfile];
D.matfileinv3x4=['3x4_' D.matfileinv];
D.origvol2almatfn=matI; % the default; could also be a file name for alignment from original anatomical to anatomical aligned to epi
D.surfvolfn=sprintf('%s_SurfVol+orig',C.sumasubjid);
D.surfvolalfn=sprintf('%s_SurfVol_al+orig',C.sumasubjid);
D.surfvol2alorigvolfn3x4='3x4_SurfVol2OrigVol.aff12.1D';
D.surfvol2alorigvolfn3x4lpirai='3x4_SurfVol2OrigVol_lpirai.aff12.1D';
D.surfalpostfix='_al';
D.afnidir='';
D.surfpat=[C.icopat '%sh.%s%s.asc']; % from afni's mapicosahedron; last '%s' is infix
D.costfunction='lpa';
D.skullstriporigvol=true; % skull strip original volume before alignment?
D.skullstripsurfvol=true; % skull strip surf volume (form freesurfer) before alignment?
D.skullstripsuffix='_ss'; % suffix for skull stripped versions of anatomical volumes
D.alignepianatsuffix='_al2SV';
D.overwrite=false;

[p,n,v]=afni_fileparts([C.origvolfn]);
if isempty(p), error('not supported - origvolfn should contain a path'); end
D.origvolsvgridfn=[n '_SVgrid' v];

C=surfing_struct(D,varargin{:});



cmd=C.shellinit;
afnidir=C.afnidir; % used by cmdaddafni (defined at end of this file)


allactions={'tosuma','mapico','moresurfs','alignoriganat','applyoriganat'};

if isnumeric(C.action)
    C.action={allactions{C.action}};
end

if iscell(C.action)
    actions=C.action;
    nactions=numel(actions);
    cmds=cell(1,nactions);
    for k=1:nactions
        C.action=actions{k};
        cmds{k}=me(C);
    end
    cmd=[cmds{:}];
else

    switch C.action
        case 'tosuma'
            cmdaddq('cd %s',C.fssurfdir);
            cmdadd('rm -rf %s',C.sumadir);
            cmdaddafni('\@SUMA_Make_Spec_FS', '-sid %s', C.sumasubjid);

        case 'mapico'
            cmdaddq('cd %s',C.sumadir)

            for hemi=C.hemis
                cmdaddafni('MapIcosahedron', ['-overwrite -spec %s_%sh.spec -ld %d -prefix ' C.icopat],...
                    C.sumasubjid,hemi,C.icold,C.icold);
            end

        case 'moresurfs'
            % average surfaces to get intermediate and semi-inflated surfaces
            n=size(C.moresurfnames,1);

            for hemi=C.hemis
                for k=1:n
                    fn1=sprintf([C.sumadir C.surfpat],C.icold,hemi,C.moresurfnames{k,1},'');
                    fn2=sprintf([C.sumadir C.surfpat],C.icold,hemi,C.moresurfnames{k,2},'');
                    fn3=sprintf([C.sumadir C.surfpat],C.icold,hemi,C.moresurfnames{k,3},'');

                    [v1,f1]=freesurfer_asc_load(fn1);
                    [v2,f2]=freesurfer_asc_load(fn2);
                    if ~isequal(f1,f2)
                        error('Non-matching topollogy for %s and %s',fn1,fn2);
                    end

                    v3=(v1+v2)/2;

                    freesurfer_asc_save(fn3,v3,f1);
                    fprintf('Averaged %s and %s -> %s (%sh)\n',C.moresurfnames{k,1},C.moresurfnames{k,2},C.moresurfnames{k,3},hemi);
                end
            end

        case 'alignoriganat'
            % find transformation through anatomical
            [p,n,v,e,e2]=afni_fileparts([C.sumadir C.surfvolfn]);
            if isempty(p)
                error('Not supported: relative path');
            end

            cmdaddq('mkdir %s; cd %s',C.aligndir,C.aligndir);
            if ~exist([C.aligndir n v e],'file') && ~exist([C.aligndir n v e2],'file') % check that surfvolfn does not exist yet
                cmdaddq('cp %s/%s%s%s %s/%s%s%s',p,n,v,e ,C.aligndir,n,v,e ); % copy head and brik
                cmdaddq('cp %s/%s%s%s %s/%s%s%s',p,n,v,e2,C.aligndir,n,v,e2);
            else
                warning('Not overwriting: %s', [p '/' n v e]);
            end

            % skull strip original and surf volumes, if necessary
            origsurfvols={[C.sumadir C.surfvolfn],C.origvolfn};
            striporigsurfvols=[C.skullstriporigvol C.skullstripsurfvol];
            alignvols=cell(1,2); % names of files to be used for alignment
            
            skullyesopts={'-anat_has_skull yes',''};
            skullnoopts={'-epi_strip None','-anat_has_skull no'};
            
            alignopts='';
            
            for k=1:2
                [p,n,v,e,e2]=afni_fileparts(origsurfvols{k});
                if striporigsurfvols(k)
                    newfn=['./' n C.skullstripsuffix v];
                    if exist([C.aligndir newfn e],'file')
                        warning('For skull strip, not overwriting %s',newfn);
                    else
                        cmdaddafni('3dSkullStrip','-overwrite -input %s -prefix %s',origsurfvols{k},newfn)
                    end
                    alignvols{k}=newfn;
                    alignopts=[alignopts ' ' skullnoopts{k}];
                else
                    alignvols{k}=[n v];
                    alignopts=[alignopts ' ' skullyesopts{k}];
                end
            end
                   
            % clean up add bit
            cmdadd('rm %s %s %s %s',C.matfileinv,C.matfile,C.matfile3x4,C.matfileinv3x4);

            % find transformation to alineated file
%            cmdaddafni('3dAllineate','-overwrite -base %s -source %s -1Dmatrix_save %s -cmass -warp shr -cost %s',C.origvolfn,C.surfvolfn,C.matfileinv,C.costfunction);
            %cmdaddafni('3dAllineate','-overwrite -automask -source_automask -base %s -source %s -1Dmatrix_save %s -cmass -warp shr -cost %s',alignvols{1},alignvols{2},C.matfileinv,C.costfunction);
            cmdaddafni('align_epi_anat.py','-overwrite -dset1 %s -dset2 %s -dset1to2 -giant_move -suffix %s -Allineate_opts "-warp shr -VERB -weight_frac 1.0" %s',...
                alignvols{2},alignvols{1},C.alignepianatsuffix,alignopts);
            
            [dummy,no]=afni_fileparts(alignvols{2});
            alignmatfn=[no C.alignepianatsuffix '_mat.aff12.1D'];
            cmdaddafni('cat_matvec','-ONELINE %s -I > %s',alignmatfn,C.matfileinv);
            %cmdaddp('cat_matvec','-ONELINE %s -I > %s',matfileinv,matfile); 

        case 'applyoriganat'
            % apply found alignment parameters to volume and surfaces
            cmdaddq('cd %s',C.aligndir);
            
            %fns={C.surfvolfn,C.surfvolfixfn;C.origvolalfn,C.origvolalfixfn};
            fns={C.surfvolfn,C.surfvolalfn};
            if C.overwrite || ~exist(fullfile(C.aligndir,[C.surfvolalfn '.HEAD']),'file')
                for k=1:size(fns,1);
                    [dummy,ns,es]=afni_fileparts(fns{k,1});
                    src=[ns es];

                    [dummy,nt,et]=afni_fileparts(fns{k,2});
                    trg=[nt et];

                    %cmdaddafni('3dcopy', '-overwrite %s %s',src,trg); 
                    % add -I %here
                    cmdaddafni('cat_matvec','-ONELINE %s %s > %s',C.origvol2almatfn,C.matfileinv,C.surfvol2alorigvolfn3x4);

                    cmdadd('dims=`3dAttribute DELTA %s | sed s/-//g`; dim1=`echo $dims | cut -f1 -d" "`; if [ "$dim1 $dim1 $dim1 " = "$dims" ]; then alopt="-newgrid $dim1"; echo "New grid with dimension $dim1 mm isotropic, using option $alopt"; else alopt=""; echo "No new grid: different dimensions $dims"; fi', src);


                    if k==1 % surfvol
    %                    cmdaddafni('3dAllineate','-overwrite -1Dmatrix_apply %s -prefix %s -source %s $alopt',C.surfvol2alorigvolfn3x4,trg,src);
                        cmdaddafni('3dWarp', '-overwrite -matvec_out2in `cat_matvec -MATRIX %s` -prefix %s $alopt %s',C.surfvol2alorigvolfn3x4,trg,src);

                    else %  - obsolete
                        cmdaddafni('3dcopy','-overwrite %s %s',src,trg);
                    end
                    headers={'ALLINEATE_MATVEC_B2S_000000','ALLINEATE_MATVEC_S2B_000000','WARPDRIVE_MATVEC_FOR_000000','WARPDRIVE_MATVEC_INV_000000'};
                    for m=1:numel(headers)
                        cmdaddafni('3drefit','-atrfloat %s %s %s',headers{m},'"1 0 0 0 0 1 0 0 0 0 1 0"',trg);
                    end
                end

                % resample the original anat_al to the same grid as surfvol
                cmdaddafni('3dresample','-overwrite -prefix ./%s -master %s -inset %s',...
                    C.origvolsvgridfn,C.surfvolalfn,C.origvolfn);

                cmdadd('rm *_e3+orig.????* *_ec+orig.????* _ae.* *_e+orig.????*');
                cmdaddafni('\@Addedge','%s %s',C.origvolsvgridfn,C.surfvolalfn);
            end
            
            cmdaddafni('cat_matvec','-ONELINE %s %s %s > %s',matlpirai,C.surfvol2alorigvolfn3x4,matlpirai,C.surfvol2alorigvolfn3x4lpirai)

            for hemi=C.hemis
                n=numel(C.allsurfnames);
                for k=1:n
                    fn1=sprintf([C.sumadir C.surfpat],C.icold,hemi,C.allsurfnames{k},'');
                    fn2=sprintf([C.aligndir C.surfpat],C.icold,hemi,C.allsurfnames{k},'');
                    fn3=sprintf([C.aligndir C.surfpat],C.icold,hemi,C.allsurfnames{k},C.surfalpostfix);
                    if C.overwrite || ~exist(fn3)
                        cmdaddq('cp %s %s',fn1,fn2);
                        cmdaddafni('ConvertSurface', '-overwrite -i_fs %s -o_fs %s -ixmat_1D %s',fn2,fn3,C.surfvol2alorigvolfn3x4lpirai);
                    end
                end
            end
            
            
        case 'all'
            C.action=allactions;
            cmd=me(C); 
            return
            
        case '',

        otherwise
            error('Did not recognize action: ''%s''', C.action);

    end
    runcmd(cmd);
end

function runcmd(cmd)

fprintf('Running: %s\n', cmd);
[s,w]=unix(cmd,'-echo');

if s~=0
    w
    error('Error running command');
end


% add to the variable 'cmd' in the caller. Syntax is like sprintf
function cmdadd(pat,varargin)
t=sprintf(pat,varargin{:});
e=sprintf('cmd=[cmd ''%s;''];',t);
evalin('caller',e);

% as above, but quit if the command fails (when executed)
function cmdaddq(pat,varargin)
t=sprintf(pat,varargin{:});
e=sprintf('cmd=[cmd ''%s||exit 1;''];',t);
evalin('caller',e);

% as above, but use the AFNI path if it is specified
function cmd=cmdaddafni(prog,pat,varargin)
afnidir=evalin('caller','afnidir');
cmd=sprintf('%s %s',fullfile(afnidir,prog),sprintf(pat,varargin{:}));
e=sprintf('cmd=[cmd ''%s||exit 1;''];',cmd);
evalin('caller',e);