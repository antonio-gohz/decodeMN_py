function [xlims, ylims] = ginput_zoom_limits(varargin)
% zoomCursor zooms in/out of an axes object and focuses on the current
% mouse position.
% zoomCursor(factor)  - Zooms by specified factor.
%                       When factor > 1, zoom in by factor.
%                       When 0 < factor < 1, zoom out by 1/factor.
% zoomCursor(fig, __) - takes the active axes object from the figure 
%                       specified by fig, and uses zoomCursor on those axes
% zoomCursor(ax, __)  - uses zoomCursor on the axes specified by ax
%
zoomAxis = 'on';

if nargin < 1  % no argument
    ax = gca();
elseif ishandle(varargin{1}) && isa(varargin{1}, 'matlab.ui.Figure')
    ax = gca(varargin{1});
elseif isa(varargin{1}, 'matlab.graphics.axis.Axes')
    ax = varargin{1};
else  % only factor as argument
    ax = gca();
end

if nargin > 1  % 
    zoomAxis = varargin(2);
    if contains(zoomAxis,'x')
        zoomAxis = 'xon';
    elseif contains(zoomAxis,'y')
        zoomAxis = 'Yon';
    end
end
zoom(ax,zoomAxis)
%waitforbuttonpress
pause;
xlims= get(ax,'Xlim');
ylims= get(ax,'Ylim');



