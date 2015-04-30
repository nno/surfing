% latex-prep: matlab script to be run before generating documentation

surfingdir=[fileparts(which('surfing_voxelselection')) '/../'];
subdirs={'surfing','misc','afni'};
trgdir='./helpfiles/';

if ~exist(trgdir,'file')
    mkdir(trgdir);
end

for j=1:numel(subdirs);
    d=[surfingdir '/' subdirs{j}];
    [fns,n]=surfing_dir([d '/*.m']);

    for k=1:n
        [p,nm,ex]=fileparts(fns{k});
        r=help(nm);
        r=regexprep(r,'%','%%');
        fnout=sprintf('%s/%s.txt',trgdir,nm);
        fid=fopen(fnout,'w');
        fprintf(fid,'%s',r);
        fclose(fid);

    end
    d
end

% symbolic link to surfing (unix-only)
if ~exist('./surfing','file')
    unix(sprintf('ln -s %s surfing',surfingdir));
end
