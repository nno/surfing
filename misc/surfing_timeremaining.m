function [x,y]=surfing_timeremaining(progress,prefix,postfix,prev_progress_line)
% Prints a status message of the time elapsed and remaining for set of
% operations/
%
% I) Old usage
%
% SURFING_TIMEREMAINING(PROGRESS,PREFIX,POSTFIX)
% INPUTS: 
%   PROGRESS    Scalar between 0 and 1 indicating how many of the 
%               operations have been completed (0,1 means none,all)
%   PREFIX   }  Text added before or 
%   POSTFIX  }             after the status message
%
% This function prints the text "took X sec, ETA Y sec", surrounded by
% prefix and postfix if they are given (ETA=Estimated Time of Arrival). If
% prefix is omitted, it is set to the calling function
%
% It is assumed that when the operations were started, the function TIC() 
% has been evoked. The present function calls TOC() to estimate the time 
% remaining using linear extrapolation.
%
% If called with one output argument, S=SURFING_TIMEREMAINING(...), then no
% message is printed but the message is returned in the string S.
%
% If called with two output arguments, [X,Y]=SURFING_TIMEREMAINING(...),
% then the number of seconds it took (X) and the ETA (Y) are returned.
%
% **************************************************************
%
% II) New usage (from CoSMoMVPA, cosmo_show_progress; added Apr 2014)
%
% progress_line=surfing_timeremaining(clock_start, progress[, msg[, prev_progress_line]])
%
% Inputs:
%   clock_start         The time the task started (from clock()).
%   progress            0 <= progress <= 1, where 0 means nothing 
%                       completed and 1 means fully completed.
%   msg                 String with a message to be shown next to the 
%                       progress bar (optional).
%   prev_progress_line  The output from the previous call to this
%                       function, if applicable (optional). If provided
%                       then invoking this function prefixes the output
%                       with numel(prev_progress_msg) backspace characters,
%                       which deletes the output from the previous call 
%                       from the console. In other words, this allows for 
%                       showing a progress message at a fixed location in 
%                       the console window.
%
% Output: 
%   progress_line       String indicating current progress using a bar,
%                       with time elapsed and expected time to complete 
%                       (using linear extrapolation). 
%
% Notes:
%   - As a side effect of this function, progress_msg is written to standard
%     out (the console). 
%   - The use of prev_progress_line may not work properly if output is
%     written to standard out without using this function.
%
% Example:
%   % this code takes just over 3 seconds to run, and fills a progress bar.
%   prev_msg=''; 
%   clock_start=clock(); 
%   for k=0:100
%       pause(.03); 
%       status=sprintf('done %.1f%%', k);
%       prev_msg=surfing_timeremaining(clock_start,k/100,status,prev_msg); 
%   end
%   % output:
%   > +00:00:03 [####################] -00:00:00  done 100.0%
%
% See also TIC, TOC
%
% NNO Dec 2010

old_usage=nargin<2 || ischar(prefix);

if old_usage
    if nargin<2, 
        [st,i]=dbstack();
        if numel(st)>=i+1
            prefix=[st(i+1).name ': ']; % caller
        else
            prefix='';
        end
    end
    if nargin<3, postfix=''; end

    tc=toc();
    eta=(1-progress)/progress*tc;

    msg=sprintf('%stook %d sec, ETA %d sec%s', prefix,round(tc), round(eta),postfix);

    switch nargout
        case 0
            fprintf('%s\n', msg);
        case 1
            x=msg;
        case 2
            x=tc;
            y=eta;
    end
else
    % translate variable names from old to new usage
    clock_start=progress;
    progress=prefix;
    msg=postfix;
    
    if nargin<4 || isempty(prev_progress_line)
        delete_count=0; % nothing to delete
    elseif ischar(prev_progress_line)
        delete_count=numel(prev_progress_line); % count the characters
    end % if not a string, die ungracefully
    
    if nargin<3 || isempty(msg)
        msg='';
    end
    if progress<0 || progress>1
        error('illegal progress %d: should be between 0 and 1', progress);
    end
    
    took=etime(clock, clock_start);
    eta=(1-progress)/progress*took; % 'estimated time of arrival'
   
    % set number of backspace characters
    delete_str=repmat('\b',1,delete_count);
    
    % define the bar
    bar_width=20;
    bar_done=round(progress*bar_width);
    bar_eta=bar_width-bar_done;
    bar_str=[repmat('#',1,bar_done) repmat('-',1,bar_eta)];
    
    % because msg may contain the '%' character (which is not to be 
    % interpolated) care is needed to ensure that neither building the 
    % progress line nor printing it to standard out applies interpolation. 
    progress_line=[sprintf('+%s [%s] -%s  ', secs2str(took), bar_str, ...
                                        secs2str(-eta)),...
                   msg,...           % avoid pattern replacement of '%'
                   sprintf('\n')];
                                      
    % the '%%' occurences in msg will be replaced back to '%' by fprintf
    fprintf([delete_str strrep(progress_line,'%','%%')]);
    x=progress_line;
    y=[];
end
    






function [m,d]=moddiv(x,y)
    % helper function that does mod and div together so that m+d*y==x
    m=mod(x,y);
    d=(x-m)/y; 

function str=secs2str(secs)
    % helper function that formats the number of seconds as 
    % human-readable string

    % make secs positive (calling function should add '+' or '-')
    secs=abs(secs);

    if ~isfinite(secs)
        str='oo'; % attempt to look like 'infinity' symbol
        return
    end

    secs=round(secs); % do not provide sub-second precision

    % compute number of seconds, minutes, hours, and days
    [s,secs]=moddiv(secs,60);
    [m,secs]=moddiv(secs,60);
    [h,d]=moddiv(secs,24);

    % add prefix for day, if secs represents at least one day
    if d>0
        daypf=sprintf('%dd+',d);
    else
        daypf=''; 
    end

    str=sprintf('%s%02d:%02d:%02d', daypf, h, m, s);
