Surfing
=======
Surfing is a Matlab toolbox for surface-based voxel neighborhood selection on the cerebral cortex, intended for informationg mapping of functional magnetic resonance imaging (fMRI) data. 

Features
--------
- voxel selection based on cortical surface reconstruction.
- support for Euclidian, Dijkstra and geodesic distance metrics.
- selection of voxels across the brain based on either a fixed radius or a fixed number of voxels.
- support for twin surfaces (FreeSurfer; pial and white) and single surfaces (Caret and BrainVoyager).

Non-features
------------
- no GUI: commands are entered in the Matlab command window, or in a matlab script.
- no full processing pipeline: preprocessing, surface reconstruction, classification, statistical analyses, and visualization are not part of the toolbox. Many other programs can be used for this, for example FSL, AFNI, SPM for preprocessing and regression; Freesurfer and Caret for surface reconstruction; several online toolboxes for classification and statistical analyses; and AFNI, Freesurfer and Caret for visualization. 

Requirements
------------
- A working installation of Matlab
- For geodesic distance metrics a working mex compiler environment may be required
- GIFTI surfaces require the GIFTI matlab toolbox: http://www.artefact.tk/software/matlab/gifti
- BrainVoyager surfaces require Neuroelf: neuroelf.net
- AFNI NIML I/O requires the AFNI matlab toolbox: http://afni.nimh.nih.gov/afni/matlab/

Developers
----------
Contact:

- Nikolaas N. Oosterhof [NNO]   <nikolaas.oosterhof at unitn.it>
- Tobias Wiestler [TW]          <tobias.wiestler at googlemail.com>
- Joern Diedrichsen [JD]        <j.diedrichsen at ucl.ac.uk>

Quick start
-----------
- add the "surfing", "misc" and "toolbox_fast_marching" directories (and their subdirectories) to the Matlab path.
- run "compile_mex" in the "surfing" root directory.
- view the surfing_voxelselection.m function, or consider the examples in the "examples" directory.

More information
----------------
see http://sourceforge.net/surfing

License
-------
Copyright (c) 2009-2014 Nikolaas N. Oosterhof, Tobias Wiestler, Joern Diedrichsen.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted 
provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of 
conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of 
conditions and the following disclaimer in the documentation and/or other materials provided with 
the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY 
WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Citations
---------
If you use this toolbox for a scientific publication, please cite:
Oosterhof, N.N., Wiestler, T, Downing, P.E., & Diedrichsen, J. (in press). A comparison of volume-
based and surface-based information mapping. Neuroimage. Available online 
http://dx.doi.org/10.1016/j.neuroimage.2010.04.270

If you use geodesic distances, please also cite:
Gabriel Peyre (2008), Toolbox Fast Marching. https://www.ceremade.dauphine.fr/~peyre/matlab/fast-marching/content.html

