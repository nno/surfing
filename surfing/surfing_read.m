function [v,f]=surfing_read(fn)
% general data read function
%
% [v[,f]]=surfing_read(fn)
%
% Inputs:
%   fn            filename of either:
%                 - ASCII surface (.asc)
%                 - GIFTI surface (.gii); requires gifti toolbox
%                 - BrainVoyager surface (.srf); requires neuroelf toolbox
%                 - NeuroImaging Markup Language dataset (.niml.dset);
%                   requires AFNI Matlab toolbox
%                 - text file data (.1D or .txt)
%
% Returns:
%   v             Px3 surface vertex coordinates, LPI (for surfaces); or
%                 PxQ data (for non-surfaces)
%   f             Rx3 face indices (for surfaces); not set for
%                 non-surfaces.
%
% Dependencies:
%  - for BrainVoyager: neuroelf, see neuroelf.net
%  - for FreeSurfer binary: freesurfer matlab library, see freesurfer.net
%  - for GIFTI: GIFTI library, http://www.artefact.tk/software/matlab/gifti
%  - for AFNI (NIML): http://afni.nimh.nih.gov/afni/matlab/
%
% Notes:
%  - the GIFTI library currently supports anatomical files only, not
%    functional data
%  - Coordinates in BV surfaces are transformed from ASR to LPI, and
%    the handedness of the faces is inverted. This way they conform better
%    to what all other packages use.
%
% Examples:
%  >> [v,f]=surfing_read('pial.asc');      % ASCII (AFNI SUMA)
%  >> [v,f]=surfing_read('pial.gii');      % GIFTI
%  >> [v,f]=surfing_read('pial.srf');      % BrainVoyager
%  >> niml=surfing_read('data.niml.dset'); % AFNI NIML
%  >> data=surfing_read('data.1D');        % } ASCII text
%  >> data=surfing_read('data.txt');       % } data
%
% See also: surfing_write
%
% NNO Feb 2014

ends_with = @(x,e) ischar(x)&&numel(x)>=numel(e)&&...
                strcmp(x(end+((1-numel(e)):0)),e);


if ends_with(fn,'.asc')
    [v,f]=freesurfer_asc_load(fn);

elseif ends_with(fn,'.gii')
    % requires gifti
    w=which('gifti');
    if isempty(w)
        error(['GIFTI is required: see '...
                'http://www.artefact.tk/software/matlab/gifti/']);
    end
    g=gifti(fn);
    v=g.vertices;
    f=g.faces;
elseif ends_with(fn,'.srf')
    % requires neuroelf
    w=which('xff');
    if isempty(w)
        error('Neuroelf is required: see neuroelf.net');
    end

    x=xff(fn);

    % BV does a lot of things completely different than other packages...
    %
    % transformation for ASR to LPI

    asr2lpi=[0 -1 0;0 0 -1;-1 0 0];

    % apply transformation
    v=x.VertexCoordinate*asr2lpi+128;

    % handedness is different 'for historical reasons' BV's authors
    % decided to do the opposite of everybody else. However because an odd
    % number of sign flips occurs in asr2lpi the orientation of the faces
    % does not have to be changed
    f=x.TriangleVertex+0; % make a copy

elseif ends_with(fn,'.niml.dset')
    w=which('afni_niml_readsimple');
    if isempty(w)
        error('AFNI matlab library is required: see neuroelf.net');
    end

    v=afni_niml_readsimple(fn);

elseif ends_with(fn,'.1D') || ends_with(fn,'.txt')
    v=textread(fn,'','commentstyle','shell');

else
    w=which('freesurfer_read_surf');
    if isempty(w)
        error('FreeSurfer Matlab is required; see freesurfer.org');
    end

    try
        [v,f]=freesurfer_read_surf(fn);
    catch exception
        error(['Unknown extension in %s, or unable to read FreeSurfer'...
                    ' binary file'], fn);
    end
end

