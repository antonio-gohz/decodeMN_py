function ax = plotMUs(signal,fs,varargin)
%PLOTMUS Summary of this function goes here
%   Detailed explanation goes here
% todo check sizes and transpose if necesary
numMUs = min(size(signal.spikeTrains));

% default params
refSignal = []; % torque, activation, force... 
fs_refSignal = 100;
tVec = []; % must be the same length as spiketrains
flagPlotSpikeTrains = false;
flagPlotDR = true;
titles = {'MN spike trains', 'MN discharge rates'};
labels = {'Time (s)', 'Motor neuron (#)', 'Mean discharge rate (pps)', 'Torque (Nm)'};
lineWidths = [0.5, 0.5, 3];
lineSpecs =  {'|','.','-'};
markerSizes = [10,10,0];
% plot colors: Bright qualitative is a colorblind friendly palette
colors = {'#4477AA'; '#EE6677'; '#228833'; '#CCBB44'; '#66CCEE'; '#AA3377'; '#BBBBBB'}; % Bright qualitative
colors = repmat(colors,ceil(numMUs/size(colors,1)),1);
colors = colors(1:numMUs);
dischargeRates =mean(signal.dischargeRates,2,'omitnan');
mnIndexTicks = {{num2str((1:numMUs)','%.0f')},{num2str(dischargeRates(:),'%.1f')}};


%% Parse optional input parameters
if ~isempty(varargin)
    if rem(length(varargin), 2) ~= 0
        error('Optional input arguments should be provided as parameter-value pairs.')
    end
    for i = 1:2:length(varargin)
        param = varargin{i};
        value = varargin{i+1};
        switch param
            case 'tVec'
                tVec = value;
            case 'dischargeRates'
                dischargeRates = value;
            case 'fs_refSignal'                
                fs_refSignal = value;
            case 'refSignal'
                refSignal = value;
            case 'flagPlotSpikeTrains'
                flagPlotSpikeTrains = value;
            case 'flagPlotDR'
                flagPlotDR = value;
            case 'titles'  
                titles = value;
            case 'labels' 
                labels = value;
            case 'lineWidths'
                lineWidths = value;
            case 'Colors'
                colors = value;            
            case 'mnIndexTicks'
                mnIndexTicks= value;
            otherwise
                error('Invalid optional input parameter: %s', param)
        end
    end
end

flagRS=isempty(refSignal);

% Plot handles
if flagPlotSpikeTrains && flagPlotDR % if plot both
    ax(1) = subplot(1, 2, 1);
    ax(2) = subplot(1, 2, 2);
else
    ax(1) = axes;
end
titles = titles(logical([flagPlotSpikeTrains,flagPlotDR]));
labels = labels(logical([1,flagPlotSpikeTrains,flagPlotDR,~flagRS]));
spikes = nan(size(signal.spikeTrains));
spikes(logical(signal.spikeTrains)) = 0.5;
plotSignals = {(0.8*spikes'+(1:numMUs)-0.4),signal.dischargeRates'./ 40 + (1:numMUs) - 0.5};
plotSignals = plotSignals(logical([flagPlotSpikeTrains,flagPlotDR]));
lineSpecs = lineSpecs(logical([flagPlotSpikeTrains,flagPlotDR,~flagRS]));
mnIndexTicks = mnIndexTicks(logical([flagPlotSpikeTrains,flagPlotDR]));

% make tVector 
if isempty(tVec)
    len = length(signal.spikeTrains);
    tVec = linspace(0,len/fs,len);
end

for i=1:length(ax)
    h = plot(ax(i),tVec,plotSignals{i},lineSpecs{i},'MarkerSize', 10);
    set(h,{'Color'}, colors)

    title(ax(i),titles{i});
    xlabel(ax(i),labels{1});
    ylabel(ax(i),labels{i+1});
    xlim(ax(i),[-0.5, max(tVec) + 0.5]);
    ylim(ax(i),[0, numMUs + 0.5]);    title(titles{i});
    yticks(ax(i),0:1:numMUs+1)

    if ~flagRS
        yyaxis(ax(i),'right')
        ax(i).YAxis(2).Color = 'k';
        plot(ax(i),tVec, refSignal, 'Color', [0.7, 0.7, 0.7, 0.5], 'LineWidth', lineWidths(end));
        ylabel(labels{end});
    end
    yticks(1:numMUs)
    yticklabels(mnIndexTicks{i});

end
linkaxes(ax,'x');


end

