function varargout = zoomginput(varargin)
% zoomginput activates ginput and allows zooming in/out of the active axis
% using the mouse wheel.
% scroll up to zoom in, scroll down to zoom out.
% all other functionality is identical to ginput.
%
% Author: TADA 2019
% scrollListener = addlistener(gcf(), 'WindowScrollWheel', @onScroll);
var = varargin;
idaxes = find(strcmp(var,'x') | strcmp(var,'y') );
if ~isempty(idaxes)
    selectAxes = var{idaxes};
    var(idaxes)=[];
else
    selectAxes='xy';
end

scrollListener = addlistener(gcf(), 'WindowScrollWheel',...
    @(src, edata) onScroll(src, edata, selectAxes));

try
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = ginput(var{:});
    title(['Selected channels: ',num2str(round(varargout{2}'))])
    delete(scrollListener);
catch ex
    delete(scrollListener);
    ex.rethrow();
end
end

function onScroll(src, edata,selectAxes)
% calculate zoom factor
defaultFactor = 1.1;
factor = (defaultFactor * abs(edata.VerticalScrollCount))^sign(edata.VerticalScrollCount);

zoomCursor(gca, factor,selectAxes);
end