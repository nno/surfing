function surfing_compile_mex()

% master function to compile mex functions
% including those into toolbox_fast_marching
% NNO June 2013

pwd_orig=pwd();
c=onCleanup(@()cd(pwd_orig)); % add handler to cd to original directory

me_dir=fileparts(which(mfilename()));
if ~isequal(pwd_orig,me_dir)
    error('%s must be run from its root directory %s',mfilename,me_dir);
end

surfing_infixes={'eucldist','invertmapping','uniqueidxsperrow'};
surfing_dir=fullfile(pwd_orig,'surfing');

n=numel(surfing_infixes);
for k=1:n
    infix=surfing_infixes{k};
    fn=sprintf('surfing_%s.c', infix);
    e=sprintf('cd %s; mex %s', surfing_dir, fn);
    eval(e);
end

e=sprintf('cd %s/toolbox_fast_marching; which(''compile_mex''), compile_mex;',pwd_orig);
eval(e)

