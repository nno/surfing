function [v,f]=surfing_read_surf(fn)
% read surface anatomy
%
% [v,f]=surfing_read_surf(fn)
%
% Input:
%   fn       filename of FreeSurfer ascii ('*.asc'), BrainVoyager ('.srf')
%            or FreeSurfer binary ('*').
%
% Outputs:
%   v        Px3 coordinates for P nodes
%   f        Qx3 node indices for Q faces
%
% Requires:
%  - for BrainVoyager: neuroelf, see neuroelf.net
%  - for FreeSurfer binary: freesurfer matlab library, see freesurfer.net
%  - for GIFTI: http://www.artefact.tk/software/matlab/gifti/
%
% Example:
%  >> [v,f]=surfing_read_surf('pial.asc');
%
% NNO Apr 2014

% helper function to look for file extension
endswith=@(e,x)~isempty(regexp(x,[regexptranslate('escape',e),'$']));


if endswith('.asc',fn)
    % freesurfer ascii surface
    [v,f]=freesurfer_asc_load(fn);
elseif endswith('.srf',fn)
    w=which('xff');
    if isempty(w)
        error('Neuroelf is required: see neuroelf.net');
    end

    s=xff(fn);
    v=s.VertexCoordinate(:,[3 2 1]); % swap handed-ness
    f=s.TriangleVertex;
elseif endswith('.gii',fn)
    w=which('gifti');
    if isempty(w)
        error(['GIFTI required: see '...
                'http://www.artefact.tk/software/matlab/gifti/']);
    end

    g=gifti(fn);
    v=double(g.vertices);
    f=double(g.faces);
else
    % any extension, try freesurfer
    try
        [v,f]=freesurfer_read_surf(fn);
    catch
        error('Unable to read %s', fn);
    end
end
