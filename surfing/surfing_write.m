function x=surfing_write(fn, v, f)
% general surface (data) write function
%
% surfing_write(fn, v[, f])
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
%  - for BrainVoyager: neuroelf, see neuroelf.net
%  - for FreeSurfer binary: freesurfer matlab library, see freesurfer.net
%  - for GIFTI: http://www.artefact.tk/software/matlab/gifti/
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
%  % 'v' is Px3 vertex coordinates
%  % 'f' is Qx3 node indices defining the triangles of the surface
%  % 'niml' is a struct with fields .data and .node_indices
%  % 'data' is a MxN array with data for M nodes
%  >> surfing_write('pial.asc', v, f);       % ASCII (AFNI SUMA)
%  >> surfing_write('pial.gii', v, f);       % GIFTI
%  >> surfing_write('pial.srf', v, f);       % BrainVoyager
%  >> surfing_write('data.niml.dset', niml); % AFNI NIML
%  >> surfing_write('data.1D',data);         % } ASCII text
%  >> surfing_write('data.txt',data);        % } data
%
% See also: surfing_read
%  
% NNO Feb 2014

ends_with = @(x,e) ischar(x)&&numel(x)>=numel(e)&&...
                strcmp(x(end+((1-numel(e)):0)),e);

if isnumeric(v)
    % do a few checks
    if nargin<3, 
        error('Require faces are third argument');
    end
    
    % check size of coordinates and faces
    [nv,three]=size(v);
    [nf,three_]=size(f);
    
    if three~=3 || three_~=3
        error('coordinates and faces should be Px3 resp Qx3');
    end
    
    % ensure faces have proper indices
    f_vec=f(:);
    if min(f_vec)<1 || max(f_vec)>nv || ~isequal(f_vec,round(f_vec))
        error('faces should be integers in range 1..%d',nv);
    end
end
            
if ends_with(fn,'.asc')
    freesurfer_asc_save(fn,v,f);
    
elseif ends_with(fn,'.gii')
    % requires gifti
    w=which('gifti');
    if isempty(w)
        error(['GIFTI is required: see '...
                'http://www.artefact.tk/software/matlab/gifti/']);
    end
    
    if isa(v,'gifti')
        g=v;
    else
        % build GIFTI structure
        params=struct();
        params.vertices=int32(v);
        params.faces=single(f);
        g=gifti(params);
    end
    
    save(g,fn);
    
elseif ends_with(fn,'.srf')
    % requires neuroelf
    w=which('xff');
    if isempty(w)
        error('Neuroelf is required: see neuroelf.net');
    end
    
    if isa(v,'xff')
        x=v;
        if ~isfield(x,'NrOfVertices') || ~isfield(x,'NrOfTriangles')
            error('xff structure does not contain surface');
        end
    else
        % build from scratch
        x=xff('new:srf');

        % allow for more options in BV (such as resampling meshes)
        x.ExtendedNeighbors = 1;
         
        % transformation for ASR to LPI
        asr2lpi=[0 -1 0;0 0 -1;-1 0 0];

        % BV does a lot of things completely different than other packages...
        % 
        % transformation for LPI to ASR
        x.VertexCoordinate=(v-128)*inv(asr2lpi);
        
        % handedness is different 'for historical reasons' BV's authors
        % decided to do the opposite of everybody else. However because an 
        % odd number of sign flips occurs in asr2lpi the orientation of the
        % faces does not have to be changed
        x.TriangleVertex=f;
        
        nv=size(x.VertexCoordinate,1);
        nf=size(x.TriangleVertex,1);
        x.NrOfVertices=nv;
        x.NrOfTriangles=nf;
        %set to grey
        x.VertexColor= [NaN(nv, 1), repmat([60, 60, 60], [nv,1])];
        
        %surfing_surface_nbrs returns the same neighbors as neuroelf's
        %TrianglesToNeighbors method, but in a different order. This seems
        %to matter for BrainVoyager. Therefore we use NeuroElf's method to
        %set neighbors when dealing with BrainVoyager surfaces
        nei = x.TrianglesToNeighbors;
        x.Neighbors = nei;
    end
    
    %recalculate normals
    %some pretty colors
    x.ConvexRGBA = [0.6 0.6 0.6 1];
    x.ConcaveRGBA =  [0.3 0.3 0.3 1]; 
   
    
    x.RecalcNormals();
    
    %ADD CURVATURE INFO
    smp = x.CurvatureMap;
    cases_concave = find(smp.Map.SMPData < 0);
    x.VertexColor(cases_concave, :) = repmat(x.ConcaveRGBA, [length(cases_concave), 1]);
    
    x.SaveAs(fn);
    
    %clearing objects leads to more stable behavior of
    %NeuroElf (at least this was true for its predecessor BVQXtools)
    x.ClearObject;
    
elseif ends_with(fn,'.niml.dset')
    % AFNI NIML
    w=which('afni_niml_readsimple');
    if isempty(w)
        error('AFNI matlab library is required: see neuroelf.net');
    end
    
    if nargin>2
        error('cannot have more than two arguments');
    end
    
    afni_niml_writesimple(v, fn);
    
elseif ends_with(fn,'.1D') || ends_with(fn,'.txt')
    % ASCII
    if ~isnumeric(v) 
        error('second argument must be numeric');
    end
    
    if nargin>2
        error('cannot have more than two arguments');
    end
    
    precision=7; % 7 decimals float
    ncol=size(v,2);
    
    % generate pattern for output
    pat=repmat(sprintf('%%.%df\t',precision),1,ncol);
    pat(end)=sprint('\n');
    
    % write results
    fid=fopen(fn,'w');
    fprintf(fid,pat,v');
    fclose(fid);
    
else
    error('Unknown extension in %s', fn);
end
                        
