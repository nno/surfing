function opt=surfing_afni_opts2string(varargin)
% translates options as used in AFNI from arguments to a string
%
% OPT=SURFING_AFNI_OPTS2STRING(...) takes the input arguments and converts
% them to a string as used in AFNI programs. 
% The input arguments can be of any type accepted by surfing_struct, e.g. a
% struct (where fieldnames serve as option names) or pairs of option names
% followed by their value. 
% Option names should not start with a  ('-')
%
% Numerical values are converted to strings, and boolean options are only 
% returned in OPT if they are true. 
%
% There is no support for values without option name, as is the case for
% some applications that use the last input argument as the input dataset.
%
% Fieldnames that end with a '_' are ignored.
% 
% Example: opt=surfing_afni_opts2string('mm',true,'prefix','out','I',2)
% returns the string '-mm  -prefix out -I 2 '
%
% NNO Mar 2011

C=surfing_struct(varargin);
ks=fieldnames(C); % key names
n=numel(ks);

opts=cell(1,n);

for j=1:n
    k=ks{j};
    v=C.(k);
    
    if k(end)=='_'
        continue
    end
    
    if islogical(v) && v
        w=''; % use fieldname as option
    elseif isnumeric(v)
        if round(v)==v
            w=sprintf('%.0f',v);
        else
            w=sprintf('%d',v);
        end
    elseif ischar(v)
        w=v;
    else
        error('Unrecognized type of option %s', k);
    end
    
    if k(1)=='-'
        k=k(2:end);
    end
    
    opts{j}=sprintf('-%s %s ',k,w);
end

opt=[opts{:}];

    


