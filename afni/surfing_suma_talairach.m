function [cmd,cmd2]=surfing_suma_talairach(sumadirs,trgdir,varargin)
% Transform surfaces and volumes to Talairach space and average them
%
% SURFING_SUMA_TALAIRACH(SUMADIRS,TRGDIR,...)
%
% INPUTS:
%   SUMADIRS    cell with directories created by SUMA
%   TRGDIR      directory where to store averaged results
%   ...         several names and prefixes as defined in the source code
%
% This function applies the Talairach transformation from Freesurfer's
% recon-all program to surfaces that were converted to AFNI SUMA format
% using @SUMA_Make_Spec_FS and MapIcosehedron. It is assumed that these
% files were created using SURFING_SUMA_ALIGN. It also applies the
% Talairach transformation to an anatomical volume.
%
% Subsequently, the surfaces and volumes are averaged and stored in TRGDIR.
% Spec files with the averaged surfaces and average are written there too.
% Additionally, two surface data .niml.dset files are written:
% - standard deviation of distance from average surface (node by node)
% - area of each node (a third of surrounding triangles)
%
% This function is currently *EXPERIMENTAL*
%
% See also SURFING_SUMA_ALIGN, SURFING_SUMA_MAKESPEC
%
% NNO Mar 2011


%sids=getvar('subjids');
%sumadirs=phoebe_map(@(x) getvar('sumadir',x),sids);
%trgdir=[getvar('glmdir') 'groupana'];

% set some defaults
D.icold=100;
D.icopat='ico%d_';
D.surftlrcpostfix='_tlrc';
D.anatall='bucketall_SurfVol';
D.anatavg='avg_SurfVol';
C=surfing_struct(D,varargin);

if ischar(sumadirs)
    sumadir=sumadirs;
    sumadirs=cell(1);
    sumadirs{1}=sumadir;
elseif ~iscell(sumadirs)
    error('Sumadirs should be a cell with directories');
end

curdir=pwd();
ifxs={'smoothwm','intermediate','pial','semiinflated','tqinflated','inflated'}; 
hemis='lr';
ni=numel(ifxs);      % number of surfaces
ns=numel(sumadirs);  %           subjects
nh=numel(hemis);     %           hemispheres (2, typically)
icostr=sprintf(C.icopat,C.icold); % prefix used in filenames

cmd=surfing_afni_runbinary('echo "Talairach transformation started"');

% keep track of which files were created
allsurffiles=cell(ns,ni,nh);
allanatfiles=cell(ns,1);

% convert data for subjects one by one
for k=1:ns
    sumadir=sumadirs{k};
    
    if ~isdir(sumadir), error('%s is not a directory',sumadir); end
    
    % go up two levels from SUMA dir, then find subject directory
    cd([sumadir '..' filesep() '..']);
    subjdir=pwd();
    cd(curdir);
    
    % get subject id from current directory
    slashpos=find(subjdir==filesep(),1,'last');
    if isempty(slashpos), error('Cannot find subject id'); end
    subjid=subjdir((slashpos+1):end);
    
    % ensure FreeSurfer knows where we are
    cmd=sprintf('%s;cd %s%s..;export SUBJECTS_DIR=`pwd`;cd %s;echo "Now in "`pwd`',...
        cmd,subjdir,filesep(),sumadir);
    
    % convert anatomical
    cmd=sprintf('%s;mri_convert -at ../../mri/transforms/talairach.xfm ../../mri/T1.mgz ./T1_tlrc.nii||exit 1',cmd);
    anatfn=sprintf('%s_SurfVol%s',subjid,C.surftlrcpostfix); % add '_tlrc' to filename
    cmd=sprintf('%s;3dcopy -overwrite T1_tlrc.nii %s+orig;rm %s+tlrc.????*; 3drefit -overwrite -view tlrc %s+orig',...
        cmd,anatfn,anatfn,anatfn);
    cmd=sprintf('%s;echo "Converted anatomical %s+tlrc (%s)"',cmd,anatfn,subjid);
    allanatfiles{k}=[sumadir anatfn '+tlrc'];
    
    % convert surfaces
    for h=1:nh % both hemispheres
        hemi=hemis(h);
        for j=1:ni % all infixes
            ifx=ifxs{j};
            pat=sprintf([sumadir filesep() icostr '%sh.%s.asc'],hemi,ifx);
            d=surfing_dir(pat); 
            if numel(d)~=1 % only continue if exactly one file found
                warning('No unique match for %s, skipping',pat);
                continue;
            end
            
            [dummy,src]=fileparts(d{1}); % get rid of extension
            trg=[regexprep(src,ifx,[ifx C.surftlrcpostfix])]; % add '_tlrc' to target filename
            
            % use Freesurfer mris_convert with Talairach transform
            cmd=sprintf('%s; mris_convert -t %s %s.asc %s.asc || exit 1; echo "Converted %s -> %s (%s)"',...
                cmd,subjid,src,trg,src,trg,subjid);
            allsurffiles{k,j,h}=[sumadir filesep() trg '.asc'];
        end
    end
    cmd=sprintf('%s;echo;echo "*** Completed Talairach transform for %s (%d/%d)***"; echo',cmd,subjid,k,ns);
end

cmd=sprintf('%s; mkdir %s; cd %s || exit 1',cmd,trgdir);

[a,w]=unix(cmd,'-echo');
%a=0;
if a~=0, error('Failed converting surfaces; quitting now'); end

% copy anatomical files
cmd2=surfing_afni_runbinary(sprintf('echo "Averaging anatomicals"; cd %s',trgdir));
anatlist=[];
anatfilesext=cell(1,2*ns);
hb={'HEAD','BRIK'}; % head and brik extensions
for k=1:ns
    src=allanatfiles{k};
    [dummy,n,e]=afni_fileparts(src);
    trg=[n e];    

    for j=1:numel(hb)
        % ensure volumes from different participants on same grid, by resampling to first volume
        cmd2=sprintf('%s; 3dresample -overwrite -master %s -prefix ./%s -inset %s || exit 1',...
            cmd2,allanatfiles{1},trg,src);
        anatfilesext{k*2+j-2}=[trg '.' hb{j}]; % list of anatomical files generated
    end
    anatlist=[anatlist ' ' trg]; % keep track of all anatomical names

end

% bucket and average anatomicals
cmd2=sprintf('%s; cd %s; 3dMean -overwrite -prefix %s %s; 3drefit -atrstring TEMPLATE_SPACE TLRC %s+tlrc', cmd2, trgdir,C.anatavg, anatlist, C.anatavg);


for k=1:numel(anatfilesext)
    cmd2=sprintf('%s; rm %s',cmd2,anatfilesext{k}); % remove individual volumes
end

[a,w]=unix(cmd2,'-echo');
%a=0;
if a~=0, error('Failed averaging anatomicals; quitting now'); end

% average surfaces (FIXME: no check if files exist)
for h=1:nh % for two hemispheres
    for j=1:ni % for each surface name (infxs)
        fprintf('Loading view %d / %d for %sh\n',j,ni,hemis(h));
        [v,f]=freesurfer_asc_load({allsurffiles{:,j,h}});
        n2f=surfing_nodeidxs2faceidxs(f'); % mapping of nodes to faces
        mv=mean(v,3); % average surfaces
        if j==1
            nverts=size(v,1);
            sv=zeros(nverts,ni); % standard deviation for surfaces
            sa=zeros(nverts,ni); % surface area
        elseif size(v,1)~=nverts
            error('Non-matching number of vertices - this should not happen');
        end
        
        dj=zeros(nverts,ns); % distances from mean for all subjects
        for k=1:ns
            vj=squeeze(v(:,:,k)); % coordinates of j-th participant
            d3=mv-vj; % difference in x,y,z for all nodes
            dj(:,k)=surfing_eucldist(d3',[0;0;0]); % distance from mean
            
            sa(:,j)=sa(:,j)+(1/ns)*surfing_surfacearea(vj,f,n2f);
        end
        sv(:,j)=std(dj,[],2); % compute standard deviation for this surface
        
        [dummy,fn,ext]=fileparts(allsurffiles{1,j,h});
        trgfn=[trgdir filesep() fn ext];
        freesurfer_asc_save(trgfn,mv,f); % write average surface
    end
    
    S=struct();
    S.labels=ifxs;
    prefix=sprintf('%s%s%s%sh_',trgdir,filesep(),icostr,hemis(h));
    
    % write standard deviation for all surfaces to a single file
    S.data=sv;
    afni_niml_writesimple(S,[prefix 'surface_avg_distance_std.niml.dset']);
    
    % write average surface area for all surfaces to a single file
    S.data=sa;
    afni_niml_writesimple(S,[prefix 'surface_avg_surface_area.niml.dset']);
    
end
    
surfing_suma_makespec('aligndir',trgdir,'surfvolalfn',[C.anatavg '+tlrc'],'icold',C.icold);    
