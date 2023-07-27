function [] = maximize(figureHandle)

% MAXIMIZE maximizes figure in figureHandle
% 
% author: Elliot Smith 
% Version Date: 20141020

% maximizing
screenDimensions = get(0,'Screensize');
set(figureHandle, 'Position', screenDimensions); % Maximize figure.

