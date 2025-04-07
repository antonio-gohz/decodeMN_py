function [spikeTrains, dischargeRates,mnPulses, meanDR, stdDR, recThres, refSignal, removedMUs, ax ] = ...
    getFiringProperties(mns, varargin)
% This function takes spike trains or motor neuron firing instances and
% calculates spike trains, discharge rates (continous, means, stds), 
% recruitment thresholds, resamples refSignal and plots the results.
%
% Inputs:
% - mns: Motoneuron firing events represented as EITHER: 
%       - spikeTrains: Matrix of 0s and 1s representing spike trains of 
%       motor neurons (MNs), where each MN is a columns; OR 
%       - mnPulses: Cell array where each element is a MN and it contains
%       the indices of the time instances where the motor neuron fired    
% - varargin: Variable input arguments (optional)
%     - tVec: Time vector 
%     - fs_MN: sampling frequency spikeTrains 
%     - dischargeRates: Pre-calculated discharge rates 
%     - refSignal: Reference signal 
%     - fs_refSignal: sampling frequency refSignal 
%     - flagPlotSpikeTrains: Flag to plot spike trains 
%     - flagPlotDR: Flag to plot discharge rates 
%     - outlierFlag: Flag to remove spikes with more than 3 MAD
%     - idsMeanDR: Interval to calculate DR in samples (if not specified 
%     the meanDR and stdDR are computed for the entire dataset length) 
%     - titles: for plots, 2-element cell array (1) spikeTrains and (2) 
%     dischargeRates default ->  {'MN spike trains', 'MN discharge rates'}
%     - labels: for plots, 4-element cell array (1) xlabel for all, (2)
%     and (3) ylabels for spikeTrains and dischargeRates, and (4) ylabel 
%     for refSignal. Default ->  {'Time (s)', 'Motor neuron (#)', ...
%     'Mean discharge rate (pps)', 'Force (N)'}
%     - lineWidths: 3-el vector: spikeTrains dischargeRates refSignal
%     - colors: n-element cell with HEX colors for the MNs   
%     - flagManualRemoval
%     - mnIndexTicks : vector with index (useful when plotting matched MNs 
%
% Outputs:
% - spikeTrains: Calculated spike trains
% - dischargeRates: Calculated discharge rates
% - averageDR: Average discharge rate
% - stdDR: Standard deviation of discharge rate
% - ax: Axes handles for the plots
% - refSignal: Resampled reference signal (if provided and if fs is different)
% - recruitmentThreshold: Vector of recruitment thresholds

% Input erros 
if nargin < 1
   error('Not enough input arguments'); 
end
if iscell(mns)
    mnPulses = mns;
    spikeTrains = [];
    nMNs = size(mnPulses, 2);
else
    spikeTrains = mns;
    mnPulses = cell(1, size(spikeTrains, 2));
    nMNs = size(spikeTrains, 2);
end

% Set default values
tVec = [];
dischargeRates = [];
refSignal = [];

fs_MN = 2048;
fs_refSignal = fs_MN;
outlierFlag=0;
colors = {'#4477AA'; '#EE6677'; '#228833'; '#CCBB44'; '#66CCEE'; '#AA3377'; '#BBBBBB'}; % Bright qualitative 
% colors =[0.2667    0.4667    0.6667
%     0.9333    0.4000    0.4667
%     0.1333    0.5333    0.2000
%     0.8000    0.7333    0.2667
%     0.4000    0.8000    0.9333
%     0.6667    0.2000    0.4667
%     0.7333    0.7333    0.7333];
flagPlotSpikeTrains = false;
flagPlotDR = false;
titles = {'MN spike trains', 'MN discharge rates'};
labels = {'Time (s)', 'Motor neuron (#)', 'Mean discharge rate (pps)', 'Force (N)'};
lineWidths = [0.5, 0.5, 3];
idsMeanDR = [];
removedMUs = [];
flagManualRemoval=0;
keepAllMUs=0; % this just in case you are reusing filters and wwant to see how bad the MN is 
mnIndexTicks = [];
% Process optional input parameters
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
            case 'fs_MN'
                fs_MN = value;
            case 'fs_refSignal'                
                fs_refSignal = value;
            case 'refSignal'
                refSignal = value;
            case 'outlierFlag'
                outlierFlag = value; 
            case 'idsMeanDR'
                idsMeanDR = value;
            case 'flagPlotSpikeTrains'
                flagPlotSpikeTrains = value;
            case 'flagPlotDR'
                flagPlotDR = value;
            case 'titles'   % todo delete
                titles = value;
            case 'labels'  % todo delete
                labels = value;
            case 'lineWidths'
                lineWidths = value;
            case 'Colors'
                colors = value;            
            case 'flagManualRemoval'
                flagManualRemoval = value;
            case 'mnIndexTicks'
                mnIndexTicks= value;
            case 'keepAllMUs' % previously'dontremoveanything'
                keepAllMUs= value;
            otherwise
                error('Invalid optional input parameter: %s', param)
        end
    end
end

%flags to compute dischargeRates spikeTrains mnPulses if not provided
flagDR=isempty(dischargeRates);
flagST=isempty(spikeTrains);
flagRS=isempty(refSignal);

emptyCells = cellfun('isempty',mnPulses);
flagPL= sum(~emptyCells)==0; % if all cells are empty 

% make tVector 
if isempty(tVec)
    if ~flagST
        len = length(spikeTrains);
    elseif ~flagRS && fs_refSignal==fs_MN
        len = length(refSignal);
    elseif ~flagPL  % creates a time vector with a maximum size given by the last spike of all MN pool + firstSpike (to be symetrical)
        len = max(cellfun(@max,mnPulses))+min(cellfun(@min,mnPulses));
    end
    tVec = linspace(0,len/fs_MN,len);    
else
    len = length(tVec);
end

% Initializing variables
if flagST
    spikeTrains = zeros(len,length(mnPulses));
end
if flagDR
    dischargeRates = NaN(len,length(mnPulses));
end
if flagManualRemoval
    flagPlotSpikeTrains=1;
    flagPlotDR=1;
end
% Plot handles
if flagPlotSpikeTrains && flagPlotDR % if plot both
    ax(1) = subplot(1, 2, 1);
    ax(2) = subplot(1, 2, 2);
elseif ~flagPlotSpikeTrains && ~flagPlotDR % if no plot 
    ax=[];
else
    ax(1) = axes;
    ax(2) = ax(1);
end
titles = titles(logical([flagPlotSpikeTrains,flagPlotDR]));
labels = labels(logical([1,flagPlotSpikeTrains,flagPlotDR,~flagRS]));


for mnCount = 1:nMNs
    if flagPL
        mnPulses{mnCount} = find(spikeTrains(:, mnCount)');
    end
    if flagDR
        dischargeRates(mnPulses{mnCount}(2:end), mnCount) = fs_MN ./ diff(mnPulses{mnCount});
    end
    idsNonPhysMUs = dischargeRates(:,mnCount) > 50;
    dischargeRates(idsNonPhysMUs, mnCount) = NaN;
    if flagST
        spikeTrains(mnPulses{mnCount},mnCount) = 1;
    end
    spikeTrains(idsNonPhysMUs, mnCount) = 0; % update pulses without
    mnPulses{mnCount} = find(spikeTrains(:, mnCount)'); % update pulses without
    
%     if outlierFlag
%         idsOutliersDR = isoutlier(dischargeRates(:,mnCount));
%         dischargeRates(idsOutliersDR, mnCount)=NaN;
%         spikeTrains(idsOutliersDR, mnCount)=0;
%         mnPulses{mnCount} = find(spikeTrains(:, mnCount)'); % update pulses without 
%     end

end
% Calculate average discharge rate and standard deviation
if isempty(idsMeanDR)
    idsMeanDR=1: size(dischargeRates,1);
end
% meanDR = zeros(size(idsMeanDR,1),size(dischargeRates,2));
% stdDR = zeros(size(idsMeanDR,1),size(dischargeRates,2));

%for j=1:size(idsMeanDR,1)
    meanDR= nanmean(dischargeRates(idsMeanDR,:));
    stdDR = nanstd(dischargeRates(idsMeanDR,:));
%end
if ~keepAllMUs
removedMUs= find(isnan(meanDR)| meanDR<3);
end
if outlierFlag
    % remove with few spikes lower percentile
    sumSpikes = sum(spikeTrains);
    outlierFewSpikes = isoutlier(sumSpikes,"percentiles",[5 100]);
    removedMUs = unique([removedMUs, find(outlierFewSpikes)]);
end

if ~flagManualRemoval && ~keepAllMUs
    dischargeRates(:,removedMUs)=[];
    spikeTrains(:,removedMUs)=[];
    mnPulses(removedMUs) = []; % update pulses without
    meanDR(removedMUs) = [];
    stdDR(removedMUs) = [];
end

if flagPlotSpikeTrains
    colors = repmat(colors,ceil(size(dischargeRates,2)/size(colors,1)),1);
    colors = colors((1:size(dischargeRates,2)),:);
    hSpikes = plot(ax(1),tVec,(0.8*spikeTrains+(1:size(dischargeRates,2))-0.4));
    set(hSpikes,{'Color'}, colors)
end

if flagPlotDR
    colors = repmat(colors,ceil(size(dischargeRates,2)/size(colors,1)),1);
    colors = colors((1:size(dischargeRates,2)),:);
    hDR=plot(ax(2),tVec, dischargeRates./ 40 + (1:size(dischargeRates,2)) - 0.5, '.', 'MarkerSize', 10);
    set(hDR,{'Color'}, colors)
end

if flagManualRemoval
    sgtitle("Clic on MU to remove, outlier MUs: " + num2str(removedMUs))
    [~,y]=ginput;
    removedMUs = round(unique([y', removedMUs]));
    dischargeRates(:,removedMUs)=[];
    spikeTrains(:,removedMUs)=[];
    mnPulses(removedMUs) = []; % update pulses without
    meanDR(removedMUs) = [];
    stdDR(removedMUs) = [];
    colors = repmat(colors,ceil(size(dischargeRates,2)/size(colors,1)),1);
    colors = colors((1:size(dischargeRates,2)),:);
    if ~isempty(spikeTrains)
    hDR=plot(ax(2),tVec, dischargeRates./ 40 + (1:size(dischargeRates,2)) - 0.5, '.', 'MarkerSize', 10);
    set(hDR,{'Color'}, colors)
    hSpikes = plot(ax(1),tVec,(0.8*spikeTrains+(1:size(dischargeRates,2))-0.4));
    set(hSpikes,{'Color'}, colors)
    end
end



% Resample reference signal if necessary
if ~flagRS && numel(refSignal) ~= numel(tVec)
    refSignal = resample(refSignal, numel(tVec), numel(refSignal));
end

%ax= unique(ax);
for i=1:length(unique(ax))
    title(ax(i),titles{i});
    xlabel(ax(i),labels{1});
    ylabel(ax(i),labels{i+1});
    xlim(ax(i),[-1, max(tVec) + 1]);
    ylim(ax(i),[0, size(mnPulses, 2) + 0.5]);    title(titles{i});
    yticks(ax(i),0:1:size(spikeTrains,2)+1)

    if ~flagRS
        yyaxis(ax(i),'right')
        ax(i).YAxis(2).Color = 'k';
        plot(ax(i),tVec, refSignal, 'Color', [0.7, 0.7, 0.7, 0.5], 'LineWidth', lineWidths(end));
        ylabel(labels{end});
    end
    if contains(titles{i},'rates') && isempty(mnIndexTicks)
        yticklabels({'',num2str(meanDR(:),'%.1f')});
    else
        yticklabels({'',num2str(mnIndexTicks(:),'%.0f')});
    end
        linkaxes(ax,'x');

end

% Calculate recruitment threshold
if ~flagRS
    %recThres = refSignal(cellfun(@(x) x(1), mnPulses))';
        idsRT = cellfun(@(x) x(logical([0,diff(x)>fs_MN/2])), mnPulses,'UniformOutput',false); %eliminates ids that follow each other 
                %recThres = cellfun(@(x) mean(mink(refSignal(x),10)), idsRT)';
        recThres = cellfun(@(x) mean(refSignal(x)), idsRT)';
else
    recThres = NaN;
end
