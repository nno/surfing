% master function to compile mex functions
% including those into toolbox_fast_marching
% NNO June 2013

surfing_infixes={'eucldist','invertmapping','uniqueidxsperrow'};
n=numel(surfing_infixes);
for k=1:n
    infix=surfing_infixes{k};
    fn=sprintf('surfing/surfing_%s.c', infix);
    eval(['mex ' fn]);
end

e='p=pwd(); cd toolbox_fast_marching; compile_mex; cd(p)';
eval(e)