function [v,f]=freesurfer_asc_load(varargin)
% Loads freesurfer ascii surface files
%
% [V,F]=FREESURFER_LOAD(FN) 
% INPUT:
%   FN    Filename of ascii freesurfer to load
% OUTPUT
%   V     Nx3 coordinates for N nodes
%   F     Px3 indices for nodes for P faces (base1)
%
% If multiple filenames are given, say Q, then V is a Nx3xQ array. F
% remains Px3, but it is verified that all input files have the same
% topology; if not, an error is thrown.
%
% If only one output argument is given, then S=FREESURFER_LOAD(FN) 
% returns a struct S with S.coords=V and S.faces=F.
%
% NNO June 2009, updated May 2010

if nargout==0
    return;
end

retstruct=nargout==1;

nsurfs=numel(varargin);
if nsurfs==1 && iscell(varargin{1})
    surffns=varargin{1};
    nsurfs=numel(surffns);
else
    surffns=varargin;
end

% If multiple filenames are given, then load them one by one
if nsurfs>1
    if retstruct
        v=cell(size(surffns));
        for k=1:nsurfs
            % return structs in a cell
            sk=freesurfer_asc_load(surffns{k});
            if k==1
                f=sk.faces;
            else
                fk=sk.faces; 
                d=f(:)-fk(:);
                if max(d)>0
                    error('Different topology for surface 1 and %d (%s)',k,surffns{k});
                end
            end
            v{k}=sk;
        end
    else
        for k=1:nsurfs
            [vk,fk]=freesurfer_asc_load(surffns{k});
            if k==1
                v=zeros([size(vk) nsurfs]); % allocate space for all surfaces
                f=fk;
            else
                d=f(:)-fk(:);
                if max(d)>0
                    error('Different topology for surface 1 and %d (%s)',k,varargin{k});
                end
            end
            v(:,:,k)=vk;
        end
    end
    return
end     

fn=varargin{1};

if ~exist(fn,'file')
    error('Input file %s does not exist!\n', fn);
end

fid=fopen(fn);

if fid==-1
    error('Error opening file %s\n', fn);
    return
end

% skip comment lines
while true
    pos=ftell(fid); % get current position
    line=fgetl(fid);
    if line(1) ~= '#'
        fseek(fid,pos,-1); %skip back to original position
        break;
    end
end

nvf=fscanf(fid,'%d',2);
nverts=nvf(1);
nfaces=nvf(2);

vs=fscanf(fid,'%f',4*nverts);
v4=reshape(vs,4,nverts);
v=v4(1:3,:)';

fs=fscanf(fid,'%d',4*nfaces);
f4=reshape(fs,4,nfaces);
f=f4(1:3,:)'+1; % add one to make base1

st=fclose(fid);
if st~=0
    error('Error closing file %s\n', fn);
end

if nargout==1
    S=struct();
    S.coords=v;
    S.faces=f;
    v=S;
end