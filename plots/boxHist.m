function [patchHandles] = boxHist(x,y, binWidth, hgapGrp, hgap )
%boxHist
%   This method uses patch objects to display histograms next to each boxplot.
% modified from: [Histograms] https://www.mathworks.com/matlabcentral/answers/822795-boxplot-and-histogram-in-one-plot
% REQUIREMENT:
% boxplotGroup(): https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup
% INPUTS:
% x : 1xn vector defining the x coordinate of n boxplots.
% y : mxn matrix of raw data for n boxplots
% Optional:
% binWidth = 0.4;  % histogram bin widths
% hgapGrp = .15;   % horizontal gap between pairs of boxplot/histograms (normalized)
% hgap = 0.06;     % horizontal gap between boxplot and hist (normalized)
%% Check inputs
if nargin < 5, hgap = 0.06;end
if nargin < 4, hgapGrp = .15;end
if nargin < 3, binWidth = 0.4;end
if nargin < 2, error ( ' need inputs x,y' ); end

%% Compute histogram counts & edges
hcounts = cell(size(y,2),2);
for i = 1:size(y,2)
    [hcounts{i,1}, hcounts{i,2}] = histcounts(y(:,i),'BinWidth',binWidth);
end
maxCount = max([hcounts{:,1}]);
%% Plot boxplotsGroup()
fig = figure();
ax = axes(fig);
hold(ax,'on')
% Convert y (mxn matrix) to 1xn cell array of mx1 vectors, required by boxplotWidths
yc = mat2cell(y,size(y,1),ones(1,size(y,2)));
xInterval = 1; %x-interval is always 1 with boxplot groups
normwidth = (1-hgapGrp-hgap)/2;
boxplotWidth = xInterval*normwidth;
% Define colors for each boxplot
colors = lines(size(y,2));
% Plot colored boxplots
bph = boxplotGroup(ax,yc,'BoxStyle','filled','Widths',boxplotWidth,'OutlierSize',3,'PrimaryLabels',compose('%d',x),'Colors',colors);
%set(findobj(bph.boxplotGroup,'-property','LineWidth'), 'LineWidth', 1) % increase line widths


%% Add vertical histograms (patches) with matching colors
xCoordinate = 1:size(y,2);  %x-positions is always 1:n with boxplot groups
histX0 = xCoordinate + boxplotWidth/2 + hgap;    % histogram base
maxHeight = xInterval*normwidth;       % max histogram height
patchHandles = gobjects(1,size(y,2));
for i = 1:size(y,2)
    % Normalize heights
    height = hcounts{i,1}/maxCount*maxHeight;
    % Compute x and y coordinates
    xm = [zeros(1,numel(height)); repelem(height,2,1); zeros(2,numel(height))] + histX0(i);
    yidx = [0 0 1 1 0]' + (1:numel(height));
    ym = hcounts{i,2}(yidx);
    % Plot patches
    patchHandles(i) = patch(xm(:),ym(:),colors(i,:),'EdgeColor','none','LineWidth',1,'FaceAlpha',.45);
    maxHeights(i) = max(height);
end

%% pdfs
% original: https://www.mathworks.com/matlabcentral/answers/513376-shift-horizontal-histogram-to-right#answer_422412
% modified because of boxplotGroup which makes all measurements relative to
% the groups (i.e., instead of having xIntervals of whatever amount, you
% have intervals 1,2,3,4,5 each number being a group 

% Compute the probability density estimates (pdfx) for each column of y
nSets = size(y,2); 
nPoints = 1000; 
points = linspace(min(y(:)), max(y(:)), nPoints); 
pdfx = cell2mat(arrayfun(@(col){ksdensity(y(:,col),points)'}, 1:nSets)); 
% Compute the interval between the x-values
xInterval = 1; %mean(diff(histX0)); % ie, = 5
% Normalize results so that all pdf values are 
% between 0 and 85% of the x-interval while maintaining 
% the relative heights across all pdf curves.
pdfxNorm = pdfx./max(pdfx) .* maxHeights;
%pdfxNorm = (pdfx - min(pdfx,[],'all')) * (xInterval*.85/range(pdfx(:)));
% NOTE: if you want each pdf curve's height to be independent 
% and *not* maintaining the relative hights across all pdf
% curves, use this line instead:
%   pdfxNorm = (pdfx - min(pdfx,[],1)) .* (xInterval*.85./range(pdfx,1));
% Horizontally offset the pdfxNorm values so each column
% corresponds to an associated x-value.
pdfxShift = pdfxNorm - min(pdfxNorm,[],1) + histX0;
% Compute y-values for the pdf curves so that the 
% curves do not extend beyond the range of the 
% data within each column. 
pdfYvals = repmat(points(:),1,size(pdfxShift,2)); 
pdfYvals(pdfYvals < min(y, [], 1) | pdfYvals > max(y, [], 1)) = NaN; 
% Add PDFs to plot
axis tight
ax.NextPlot = 'add';
hold on, plot(pdfxShift, pdfYvals)%'r-')
xlim([min(xCoordinate)-0.5*mean(diff(xCoordinate)),(max(histX0)+maxHeight)])


end

