function [R,S,T]=surfing_struct(varargin)
% Creates a struct from structs, pairs of keys and values, and/or cells
%
% R=surfing_struct(A1,A2,...)
%
% The "allowed form" of arguments A1,A2,... can be:
% - Ai, Ai+1 are a KEY, VALUE, pair, where KEY is a string (that can be 
%   used as a variable name) and VALUE can be of any type
% - Ai is a struct, where fieldnames are treated as KEY and the
%   corresponding value Ai.(KEY) as VALUE
% - Ai is a cell with elements in the "allowed form" (see above)
% 
% The result R is a struct with KEYs (the fieldnames) and corresponding 
% VALUEs assigned. If a KEY is used multiple times, the corresponding VALUE
% from the *last* occurence of KEY is taken.
% 
% Special modifier arguments are:
% - if Ai=='!', then Ai+1 is a cell with strings of allowed  keys
% - if Ai=='?', "                                 " required "  "
% - if Ai=='-', "                                 " to remove "  "
% Note that the use of '!', '?' and '-' are only allowed once (for now).
%
% Alternative calls:
% - without output arguments: each of the keys are assigned the value in
%   the calling workspace
% - [KEYS,VALUES,N]=surfing_struct(...) returns Nx1 cells KEYS and VALUES
%
% Examples:
%
%  R=surfing_struct('a',1,'b',true,'c','hello') returns a struct R with
%     R.a=1, R.b=true and R.c='hello'
%
%  R=surfing_struct('a',1,{'b',false,'c','hello'},b,true) does the same.
%     Note that the value of b=false is overwritten
%
%  surfing_struct('a',1,'b',true,'c','hello') assigns the values 1, true and
%     'hello' to variables a, b and c in the calling function.
%
%  S=surfing_struct(D,varargin,'?',{'a','b'}), used in a function body where D 
%     is a struct with default values, assigns the values in D first, then
%     overwrites them by the values passed to the function. If the keys 'a' 
%     and 'b' are not present (in D or varargin), an error is given.
%
% NNO Jan 2011, after processparameters (Feb 2010)

me=str2func(mfilename()); % make immune to renaming

n=numel(varargin);

Rs=cell(1,n); % cell for gathering structs (with the keys and values)
freepos=1;    % first free cell in Rs

specialchars='!?-'; % allowed and required key names, and fields to be removed
                   % note that these are not valid variable names, so their 
                   % use is not ambigous                 
specialfields=cell(1,numel(specialchars)); % support for '!', '?' and '-'                  

k=1; % argument position

while k<=n % loop over arguments, collect keys and values in Rs
    argk=varargin{k};
    if iscell(argk)
        valk=me(argk{:}); % call recursively
        
        kplus=1; % increase of k
        fplus=1; % increase of freepos
    elseif ischar(argk) && ~isempty(argk)
        if k<n
            pfind=find(argk(1)==specialchars);
            if ~isempty(pfind) % special character '?' or '!' or '-'
                if ~isempty(specialfields{pfind})
                    error('Cannot have more than one ''%s''', specialchars(pfind));
                end
                kplus=2;
                fplus=0;
                
                valk=varargin{k+1};
                if ~iscell(valk)
                    error('Value after ''%s'' should be a cell', specialchars(pfind));
                end
                specialfields{pfind}=valk;
            else % normal fieldname, make a struct with one field
                s=struct();
                s.(argk)=varargin{k+1};
                valk=s;
                kplus=2;
                fplus=1;
            end
        else
            error('Missing value after %s', argk);
        end
        
    elseif isstruct(argk) % easy - just proceed
        valk=argk;
        kplus=1;
        fplus=1;
    else
        disp(argk);
        error('The above value at argument position %d was not understood; it should be a string, struct, or cell',k);
    end
    switch fplus
        case 0, % special character, do not assign value
        case 1
            Rs{freepos}=valk; % assign value
            freepos=freepos+1; % increase freepos
        otherwise error('this should not happen');
    end
    k=k+kplus; % increase position to look for argument
end

% merge keys and values in Rs to a single struct R, in the correct order
R=struct();
for k=1:(freepos-1)
    Rsk=Rs{k};
    if isstruct(Rsk)
        fns=fieldnames(Rsk);
        for j=1:numel(fns)
            R.(fns{j})=Rsk.(fns{j});
        end
    else 
        error('This should not happen');
    end
end

minuspos=find(specialchars=='-');
toremove=specialfields{minuspos};
for k=1:numel(toremove)
    if isfield(R,toremove{k})
        R=rmfield(R,toremove{k});
    else
        warning('surfing_struct::field %s does not exist, cannot remove it',toremove{k});
    end
end

% get all keys and values in cells
keys=fieldnames(R);
n=numel(keys);
vals=cell(size(keys));
for k=1:n
    vals{k}=R.(keys{k});
end

% treat special characters '?' and '!'
for j=1:numel(specialchars)
    sfs=specialfields{j};
    if ~isempty(sfs)
        for k=1:numel(sfs)
            if ~ischar(sfs{k})
                error('Only strings in cell after ''%s'' are supported', specialchars(j)); 
            end
        end
    
        switch specialchars(j)
            case '!'
                trgs=keys;
                matches=sfs;
                msg='Not allowed';
            case '?'
                trgs=sfs;
                matches=keys;
                msg='Missing';
            case '-'
                continue
                % do nothing, should be done above
            otherwise
                error('this should not happen');
        end
           
        for k=1:numel(trgs)
            trgk=trgs{k};
            if isempty(strmatch(trgk,matches))
                error('%s field: %s', msg, trgk);
            end
        end
    end
end

% set output arguments
switch nargout
    case 0
        for k=1:n
            assignin('caller',keys{k},vals{k});
        end
    case 1
        % do nothing
    case {2,3}
        R=keys;
        S=vals;
        T=n;
    otherwise
        error('Not supported: %d output arguments', nargout);
end        