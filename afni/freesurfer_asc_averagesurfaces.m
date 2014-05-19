function freesurfer_asc_averagesurfaces(fnout,varargin)
% take the average of a number of freesurfer ASCII surface files
%
% freesurfer_asc_averagesurfaces(FNOUT,FNIN1,FNIN2,...)
% INPUTS:
%   FNOUT:    surface output file name
%   FNIN*:    surface input file name
%
% 

surffns=varargin;
surfcount=numel(surffns);

% load each surfaces
for k=1:surfcount
    surffn=surffns{k};
    if strcmpi(fnout,surffn)
        error('Output and input filename identical, will not proceed');
    end
    
    fprintf('Loading %d / %d: %s\n',k,surfcount,surffn);
    surf=freesurfer_asc_load(surffn);
    
    coords=surf.coords;
    if k==1
        % first surface, allocate space for average coordinates and store
        % faces
        avgcoords=zeros(size(coords));
        faces1=surf.faces;
    else
        % check whether this surface is same as the first one
        if ~isequal(faces1,surf.faces)
            error('Different topology for surface 1 and %d, cannot average',k);
        end
    end
    
    avgcoords=avgcoords+coords;
end

surf.coords=avgcoords / surfcount;

fprintf('Writing average of %d surfaces to %s\n',surfcount,fnout);
freesurfer_asc_save(fnout,surf);
    
