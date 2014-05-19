function [cmd,r]=surfing_afni_runbinary(cmd)
% Quick wrapper to either run a binary or get the command to run it
%
% FULLCMD=SURFING_AFNI_RUNBINARY(CMD)
%  
% INPUT:
%   CMD        Command to run
% OUTPUT:
%   FULLCMD    Command preprended with some initialization stuff; currently
%              only running .bash_profile to ensure that paths are set up
%
% If the output argument is omitted, then FULLCMD is ran with the builtin
% UNIX command.
%
% NNO Jan 2011
r=NaN;
if nargin==0, cmd=''; end

preamble='. ~/.bash_profile; . ~/.bash_login; . ~/.login; . ~/.bash_login; ~/.cshrc; ';

cmd=[preamble cmd];

if nargout==0
    s=unix(cmd,'-echo');
end