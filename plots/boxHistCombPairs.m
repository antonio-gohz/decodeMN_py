function [patchHandles] = boxHistCombPairs(x,y, binWidth, normwidthBoxPlot, hgapBoxDist, sizeHist )
%boxHistComb
%   This method uses patch objects to display histograms next to each boxplot.
% the difference between this and box Hist is that this one plots the
% distributions together
% modified from: [Histograms] https://www.mathworks.com/matlabcentral/answers/822795-boxplot-and-histogram-in-one-plot
% and for the PDF: https://www.mathworks.com/matlabcentral/answers/513376-shift-horizontal-histogram-to-right#answer_422412
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
if nargin < 6, sizeHist = 1.5;end
if nargin < 5, hgapBoxDist = 1;end
if nargin < 4, normwidthBoxPlot = 1;end
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
%% Adding a column every third column (for the distributions)
% Convert y (mxn matrix) to 1xn cell array of mx1 vectors, required by boxplotWidths
n = 3;  % for n=2, every 2nd col will be 0s
% Determine the number of columns of resultant matrix
nCol = floor(size(y,2)/(n-1)*n); 
% matrix of all NaNs
y3col = NaN(size(y,1),nCol); 
% Determine column indices of non-zeros
colIdx = 1:nCol; 
colIdx(n:n:end) = []; 
% Insert values from A into result matrix
y3col(:,colIdx) = y;

yc = mat2cell(y3col,size(y3col,1),ones(1,size(y3col,2)));
% creating extra group for distributions
xc= 1:length(yc); % NEW MOD
xInterval = 1; %x-interval is always 1 with boxplot groups

boxplotWidth = xInterval*normwidthBoxPlot;
% Define colors for each boxplot blue and red then zeros for every 3rd
% column
colors = [lines(2);0 0 0 ];
colors = repmat(colors,size(y,2)/2,1);
% Plot colored boxplots
bph = boxplotGroup(ax,yc,'BoxStyle','filled','Widths',boxplotWidth,'OutlierSize',3,'PrimaryLabels',compose('%d',xc),'Colors',colors);
%set(findobj(bph.boxplotGroup,'-property','LineWidth'), 'LineWidth', 1) % increase line widths


%% Add vertical histograms (patches) with matching colors
xCoordinate = 1:size(y,2);  %x-positions is always 1:n with boxplot groups
xGapDistributions = 3:3:size(yc,2); % gaps for plotting distributions
xGapDistributions=repelem(xGapDistributions,2);
colors = lines(2);
colors = repmat(colors,size(y,2)/2,1);
maxHeight = xInterval*normwidthBoxPlot;       % max histogram height
patchHandles = gobjects(1,size(y,2));
for i = 1:size(y,2)
    % Normalize heights
    height = hcounts{i,1}/maxCount*maxHeight * sizeHist; % MOD
    % Compute x and y coordinates
    % to plot after every pair (third column)
    xm = [zeros(1,numel(height)); repelem(height,2,1); zeros(2,numel(height))] + xGapDistributions(i)+ hgapBoxDist ; % MOD
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
pdfxShift = pdfxNorm - min(pdfxNorm,[],1) + xGapDistributions + hgapBoxDist; %MOD
% Compute y-values for the pdf curves so that the 
% curves do not extend beyond the range of the 
% data within each column. 
pdfYvals = repmat(points(:),1,size(pdfxShift,2)); 
pdfYvals(pdfYvals < min(y, [], 1) | pdfYvals > max(y, [], 1)) = NaN; 
% Add PDFs to plot
axis tight
ax.ColorOrder= colors; % to keep the same order of colors 
ax.NextPlot = 'add';  % to keep the same order of colors 
hold on, plot(pdfxShift, pdfYvals)%'r-')
%plot area under the curve 
for i = 1:size(y,2)
    patch(pdfxShift(~isnan(pdfYvals(:,i)),i), pdfYvals(~isnan(pdfYvals(:,i)),i), colors(i,:), 'FaceAlpha',0.4,'EdgeColor','none')%'r-')
end
xlim([min(xCoordinate)-0.5*mean(diff(xCoordinate)),( max(xGapDistributions)+hgapBoxDist+maxHeight*sizeHist*1.1)]) %MOD



function Jitter(data, pos, col, marker_size, line_width, cap_size)
    
    [density, value] = ksdensity(data);          
    density = density(value >= min(data) & value <= max(data));
    value = value(value >= min(data) & value <= max(data));
    value(1) = min(data);
    value(end) = max(data);

    width = 0.05/max(density);
    jitterstrength = interp1(value, density*width, data);
    jit = 2*(rand(size(data))-0.5);

    scatter(pos + jit.*jitterstrength, data, marker_size, 'filled','MarkerEdgeColor','None','MarkerFaceColor',col,'MarkerFaceAlpha',0.6);    
    errorbar(pos+0.20,mean(data),std(data),'Color','k','LineWidth',line_width,'CapSize',cap_size);
    scatter(pos+0.20, mean(data),marker_size,'filled','MarkerFaceColor','k','MarkerEdgeColor', 'None');

end 

end

