% master function to compile mex functions
% including those into toolbox_fast_marching
% NNO June 2013

pwd_=pwd();
me_dir=fileparts(which(mfilename()));
if ~isequal(pwd_,me_dir)
    error('%s must be run from its root directory %s',mfilename,me_dir);
end

surfing_infixes={'eucldist','invertmapping','uniqueidxsperrow'};
surfing_dir=fullfile(pwd_,'surfing');

n=numel(surfing_infixes);
for k=1:n
    infix=surfing_infixes{k};
    fn=sprintf('surfing_%s.c', infix);
    e=sprintf('cd %s; mex %s', surfing_dir, fn);
    eval(e);
end


e=sprintf('cd %s/toolbox_fast_marching; compile_mex;',pwd_);
eval(e)

cd(pwd_);