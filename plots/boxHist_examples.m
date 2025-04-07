%% Example boxHist
rng('default') % for reproducibility
x = 1:1:4; 
y = (randn(300, numel(x)) + linspace(.5,5,numel(x))) .* linspace(.5,2,numel(x));

binWidth = 0.4;  % histogram bin widths
hgapGrp = .05;   % horizontal gap between pairs of boxplot/histograms (normalized)
hgap = 0.2;      % horizontal gap between boxplot and hist (normalized)

[patchHandles] = boxHist(x,y, binWidth, hgapGrp, hgap );

%% Example boxHistComb
rng('default') % for reproducibility
x = 1:1:4; 
y = (randn(300, numel(x)) + linspace(.5,5,numel(x))) .* linspace(.5,2,numel(x));
binWidth = 0.4;  % histogram bin widths

sizeHist = 1.5; % gain for the size of the histogram
hgapBoxDist = 1; % distance between boxplots and distributions
normwidthBoxPlot = 1; % normalized widht of boxplot (this is related with how close the boxplots are between each other (i.e., the the bigger the widht the closer they are)
[patchHandles] = boxHistComb(x,y, binWidth, normwidthBoxPlot, hgapBoxDist, sizeHist );

%% Example boxHistCombPairs
rng('default') % for reproducibility
x = 1:1:4; 
y = (randn(300, numel(x)) + linspace(.5,5,numel(x))) .* linspace(.5,2,numel(x));
binWidth = 0.4;  % histogram bin widths

sizeHist = 1.5; % gain for the size of the histogram
hgapBoxDist = -0.5; % distance between boxplots and distributions
normwidthBoxPlot = 0.8; % normalized widht of boxplot (this is related with how close the boxplots are between each other (i.e., the the bigger the widht the closer they are)
[patchHandles] = boxHistCombPairs(x,y, binWidth, normwidthBoxPlot, hgapBoxDist, sizeHist );
legend('Pair 1', 'Pair 2')


%% Example Jitter
data = [randn(1,50) randn(1,50)+1.6 randn(1,50)+1.4  randn(1,50)+1.8]; 
cats = [cellstr(repmat('GrA',50,1)); cellstr(repmat('Group B',50,1)); cellstr(repmat('Group C',50,1)); cellstr(repmat('Group D',50,1))];   
figure
Fig = jitter_distribution_figure(data, cats, 'YLabel', 'Y label','PlotType','Internal' ); 
set(Fig.jit,'FontSize',16,'linewidth',1,'box','off')
xticks(2:3:6)
xticklabels({'Sub1','Sub2'})
xlim([0,7])

%% Example Jitter
data = [randn(1,50) randn(1,50)+1.6 randn(1,50)+1.4  randn(1,50)+1.8]; 
cats = [cellstr(repmat('GrA',50,1)); cellstr(repmat('Group B',50,1)); cellstr(repmat('Group C',50,1)); cellstr(repmat('Group D',50,1))];   
figure
Fig = jitter_distribution_pairs(data, cats, 'YLabel', 'Y label','PlotType','Internal' ); 
set(Fig.jit,'FontSize',16,'linewidth',1,'box','off')
xticks(2:3:6)
xticklabels({'Sub1','Sub2'})
xlim([0,7])