function freesurfer_asc_save(fn, v, f)
% Saves a surface to freesurfer ASCII format
%
% FREESURFER_ASC_SAVE(FN,V,F) 
% INPUTS:
%    FN:  Filename to save surface to
%    V:   3xP (or Px3) vertex coordinates
%    F:   3xQ (or Qx3) indices of faces referring to vertex coordinates
%           (F can be either base0 or base1; this is detected automatically)
%
% Alternative usage is FREESURFER_ASC_SAVE(FN,S), where S is a struct with
% fields S.coords=V and S.faces=F.
%
% NNO June 2009, updated June 2010.

if isstruct(v)
    if isfield(v, 'coords') && isfield(v, 'faces')
        verts=v.coords;
        faces=v.faces;
    else
        error('Unrecognized struct, expected fields .coords and .faces');
    end
elseif isnumeric(v) && nargin>=3 && isnumeric(f)
    verts=v;
    faces=f;
else
    error('Illegal input');
end

% transpose if necessary
[p,q]=size(verts);
if p>q, verts=verts'; end
[p,q]=size(faces);
if p>q, faces=faces'; end

[three_,nverts]=size(verts);
[three__,nfaces]=size(faces);

if three_ ~= 3 || three__ ~= 3
    error('Coordinates and faces should be Px3 and Qx3');
end

if ~isequal(faces,floor(faces))
    error('Faces should contain integers');
end

% simple check on topology
minmax=[min(faces(:)) max(faces(:))];

if isequal(minmax, [1 nverts])
    faces=faces-1; %convert base1 to base0 indexing
elseif ~isequal(minmax, [1 nverts]-1)
    error('Illegal face definition; check number of vertices\n');
end

% open file for writing
fid=fopen(fn,'w');
if fid==-1
    error('Error opening file for writing: %s\n', fn);
    return
end

% write headers
headerstr=sprintf('# Saved by freesurfer_asc_save on %s to %s\n',datestr(clock()),fn);
fprintf(fid,headerstr);
fprintf(fid,'%.0f %.0f\n', nverts, nfaces);

% write coordinates and faces
fprintf(fid,'%.6f  %.6f  %.6f  0\n', verts(:));
fprintf(fid,'%.0f %.0f %.0f 0\n', faces(:));

% close file
st=fclose(fid);
if st~=0
    error('Error closing file %s\n', fn);
end