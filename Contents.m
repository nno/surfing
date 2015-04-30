% "surfing" toolbox for surface-based voxel selection, useful for
% SURFace-based information mappING" of the cerebral cortex.
%
% Authors: Nikolaas N. oosterhof [NNO]   <n.oosterhof@bangor.ac.uk>
%          Tobias Wiestler [TW]          <tobias.wiestler.09@acl.ac.uk>
%          Joern Diedrichsen [JD]        <j.diedrichsen@acl.ac.uk>
%
% Version 0.4, May 2014
% Copyright 2009-2014
%
% If you use this toolbox for a scientific publication, please cite:
% Oosterhof, N.N., Wiestler, T, Downing, P.E., & Diedrichsen, J. (in press)
% A comparison of volume-based and surface-based information mapping.
% Neuroimage. DOI:10.1016/j.neuroimage.2010.04.270
%
% This toolbox uses the Fast Marching toolbox by Gabriel Peyre (2009),
% http://www.ceremade.dauphine.fr/~peyre/download/
%
%
% === TOOLBOX CONTENTS ===
%
% Main function
%  SURFING_VOXELSELECTION     - Voxel selection using cortical surfaces
%                               for euclidian and geodesic distances
%
% Helper functions
%  SURFING_CIRCLEROI           - Construct a circular ROI on a surface
%  SURFING_COORDS2LINVOXELIDXS - Map from coordinates to linear voxel indices
%  SURFING_EUCLDIST            - Euclidian distance between points
%  SURFING_DIJKSTRADIST        - Dijkstra distance between points
%  SURFING_NODEIDXS2COORDS     - Map from vertex indices to coordinates
%  SURFING_NODEIDXS2FACEIDXS   - Map from vertex indices to face indices
%  SURFING_READ                - Read surfaces or surface data
%  SURFING_SELECTKFIRSTIDXS    - Selects neareast K voxel indices
%  SURFING_SUBSURFACE          - Generates a 'smaller' surface based on a
%                               'bigger' surface using euclidian distance
%  SURFING_UNIQUEIDXSPERROW    - Remove duplicate indices in each row
%  SURFING_WRITE               - Write surfaces or surface data
%
% Utilities
%  SURFING_INVERTMAPPING       - inverts a matrix-stored mapping f:N->P(N)
%  SURFING_INDS2SUBS           - linear to sub indices
%  SURFING_SUBS2INDS           - sub to linear indices
%  SURFING_VOXEL2NODE          - which voxels are near which node
%
% Miscellaneous (misc directory)
%  FREESURFER_ASC_LOAD         - Reads Freesurfer ASCII surfaces
%  FREESURFER_ASC_SAVE         - Writes Freesurfer ASCII surfaces
%  SURFING_CLUSTERAREA         - Computes the area of clusters
%  SURFING_CLUSTERTFCE         - Cluster-free Threshold Enhancement (experimental)
%  SURFING_LOWRES2VONOROI      - Partition surface in patches
%  SURFING_MAPLOW2HIRES        - Correspondence between hi-low res surfaces
%  SURFING_RANDOMSUBSPACE_SELECT - Selects subsets of voxels
%  SURFING_REDUCEMAPPING       - Find which voxel data has to be loaded
%  SURFING_SPHERICAL_MASK      - Volume-based mask
%  SURFING_STRUCT              - Merge structs and input arguments
%  SURFING_SUBSAMPLE_SURFACE   - Resample surface to lower resolution
%  SURFING_SUBSPHERE2RECT      - Map part of spherical surface to rectangle
%  SURFING_SURFACE_NBRS        - Find neighbours of nodes
%  SURFING_SURFACEAREA         - Determine area of nodes and faces
%  SURFING_TIMEREMAINING       - Estimates how much time is remaining
%  SURFING_VOXELSELECTION_VOLUME - Quick volume-based voxel selection
%
% AFNI specific (afni directory)
%  AFNI_FILEPARTS              - Fileparts extended for AFNI files
%  SURFING_AFNI2SPMVOL         - AFNI info struct converted to SPM
%  SURFING_AFNI_CONJUNCTION    - Conjunction (minimum) across a dimension
%  SURFING_AFNI_GLM            - Quick GLM with t-test
%  SURFING_AFNI_OPTS2STRING    - Converts AFNI options to a string
%  SURFING_AFNI_OPTS_SUBVOL    - Correct info struct after loading sub-volume
%  SURFING_AFNI_RUNBINARY      - Runs binary after setting paths
%  SURFING_AFNI_SURFSMOOTH     - Wrapper for SurfSmooth
%  SURFING_AFNI_VOL2SURF       - Wrapper for 3dVol2Surf
%  SURFING_DIR                 - Lists files matching pattern
%  SURFING_SUMA_ALIGN          - Runs alignment for SUMA to aligned anat
%  SURFING_SUMA_MAKESPEC       - Makes SUMA .spec file
%  SURFING_SUMA_SURFACEFILES   - Lists surface and supporting files
%  SURFING_SUMA_TALAIRACH      - Talairach alignment and averaging.
