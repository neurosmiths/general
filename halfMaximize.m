function [] = halfMaximize(figureHandle,side)

% HALFMAXIMIZE maximizes figure in figureHandle to half the horizontal area
% of the currently used monitor.
%
% halfMaximize(gcf,'right')
%
% author: Elliot Smith
% Version Date: 20170530

% default side
if ~exist('side','var')
    side = 'left';
end

switch side
    case {'page'}
        set(figureHandle,'Units','Inches');
        pos = get(figureHandle,'Position');
        set(figureHandle, 'Position', [0.5 0 8.5 11]);
    case {'left'}
        % maximizing
        screenDimensions = get(0,'Screensize');
        screenDimensions(3) = ceil(screenDimensions(3)./2);
        set(figureHandle, 'Position', screenDimensions); % Maximize figure.
    case {'right'}
        % maximizing
        screenDimensions = get(0,'Screensize');
        screenDimensions(3) = ceil(screenDimensions(3)./2);
        screenDimensions(1) = floor(screenDimensions(3)./2);
        set(figureHandle, 'Position', screenDimensions); % Maximize figure.
end
