function [] = updateUser(operationMessage,loopVariable,period,NloopVariables)
% UPDATEUSER displays progress through a for loop in the command window.
%
%   [updateMessage] = updateUser(operationMessage,loopVariable,period,NloopVariables)
%
% operationMessage: a message about what kind of operation is being
%   performed (e.g. calculating spectrograms) 
% loopVariable: the loop variable
% period: number of loop variables after which to update the user
% NloopVariables: max number of loop variables

% author EHS20170609
% updated 20181010: more intuitive input args and more efficient

if isequal(mod(loopVariable,period),0)
    fprintf('\n %s %d out of %d',operationMessage,loopVariable,NloopVariables)
end

