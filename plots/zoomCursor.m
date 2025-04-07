function zoomCursor(varargin)
% zoomCursor zooms in/out of an axes object and focuses on the current
% mouse position.
% zoomCursor(factor)  - Zooms by specified factor.
%                       When factor > 1, zoom in by factor.
%                       When 0 < factor < 1, zoom out by 1/factor.
% zoomCursor(fig, __) - takes the active axes object from the figure 
%                       specified by fig, and uses zoomCursor on those axes
% zoomCursor(ax, __)  - uses zoomCursor on the axes specified by ax
%
% Author: TADA 2019
    args = varargin;
    if nargin < 1  % no argument
        ax = gca();
    elseif ishandle(varargin{1}) && isa(varargin{1}, 'matlab.ui.Figure')
        ax = gca(varargin{1});
        args = varargin(2:end);
    elseif isa(varargin{1}, 'matlab.graphics.axis.Axes')
        ax = varargin{1};
        args = varargin(2:end);
    else  % only factor as argument 
        ax = gca(); 
    end
zoomX=1;
zoomY=1;
if nargin > 2  % no argument
    zoomAxis = varargin(3);
    if strcmp(zoomAxis,'x')
        zoomY = 0;
    elseif strcmp(zoomAxis,'y')
        zoomX =0;
    end
end

    
    opt = parseZoomCursorArgs(args);
    
    axisCoord = get(ax, 'CurrentPoint');
    xyCoord = axisCoord(1, 1:2);
    if zoomX
        % calc x limits
        xlims = ax.XAxis.Limits;
        dx = xyCoord(1) - xlims(1);
        width = range(xlims);
        xlims(1) = xyCoord(1) - opt.Factor * dx;
        xlims(2) = xlims(1) + opt.Factor * width;
        ax.XAxis.Limits = xlims;
    end
    if zoomY
        % calc y limits
        for i = 1:numel(ax.YAxis)
            ylims = ax.YAxis(i).Limits;
            dy = xyCoord(2) - ylims(1);
            height = range(ylims);
            ylims(1) = xyCoord(2) - opt.Factor * dy;
            ylims(2) = ylims(1) + opt.Factor * height;
            ax.YAxis(i).Limits = ylims;
        end
    end
end

function res = parseZoomCursorArgs(args)
    res = struct();
    assert(numel(args) >= 1, 'Must specify zoom factor, by how much do you want to zoom in/out?');
    x = args{1};
    assert(isnumeric(x) && numel(x) == 1 && ~isnan(x) && isfinite(x) && x > 0 && isreal(x), 'zoomCursor factor must be a positive real finite numer. Factor value: %d', x)
    res.Factor = x;
%     parser = inputParser();
%     parser.CaseSensitive = false;
%     parser.FunctionName = 'zoomCursor';
%     parser.addRequired('Factor', 1.1,...
%         @(x) assert(isnumeric(x) && numel(x) == 1 && ~isnan(x) && isfinite(x) && x > 0 && isreal(x), 'zoomCursor Factor must be a positive real finite numer. Factor value: %d', x));
%     parser.parse(args{:});
%     res = parser.Results;
end