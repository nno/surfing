function surfing_set_path()
% helper to set the path for surfing functions

root_dir=fullfile(fileparts(which(mfilename)),'..');
sub_dirs={'misc','surfing','afni'};

for k=1:numel(sub_dirs)
    addpath(fullfile(root_dir,sub_dirs{k}));
end

fast_marching_dir=fullfile(root_dir,'toolbox_fast_marching');
addpath(genpath(fast_marching_dir));

