function plotSummaryQualityMetrics(qualitySummary, metricLabels,xLabels, varargin)

metricThresholds = nan(size(metricLabels));
vlineSeparation=[];
titlePlot =[];

% Process optional input parameters using a loop and switch case
for j = 1:2:length(varargin)
    param = varargin{j};
    value = varargin{j+1};
    switch param
        case 'metricThresholds'
            metricThresholds = value;
        case 'vectorSeparationSubjects'          
            vlineSeparation = value;
        case 'titlePlot'   
            titlePlot = value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end

subplotNumber = length(metricLabels);

figure;
x_outlier= 1:size(qualitySummary,1);
for i=1:subplotNumber
    ax(i)=subplot(subplotNumber,1,i);
    var=qualitySummary(:,i);
    [TF,L,U,C] = isoutlier(var);
    plot(x_outlier,var)
    hold on
    if contains(lower(metricLabels{i}),'sil') || contains(lower(metricLabels{i}),'pnr') ||...
        contains(lower(metricLabels{i}),'snr') || contains(lower(metricLabels{i}),'q')||...
        contains(lower(metricLabels{i}),'psd') || contains(lower(metricLabels{i}),'MUs')
        if isnan(metricThresholds(i))
            plot(x_outlier(var<L),var(var<L),'or')
        else
            plot(x_outlier(var<metricThresholds(i)),var(var<metricThresholds(i)),'or')
        end
    elseif contains(lower(metricLabels{i}),'var') || contains(lower(metricLabels{i}),'std') ||...
        contains(lower(metricLabels{i}),'cov') 
        if isnan(metricThresholds(i))
            plot(x_outlier(var>U),var(var>U),'or')
        else
            plot(x_outlier(var>metricThresholds(i)),var(var>metricThresholds(i)),'or')
        end
        plot(x_outlier(var>U),var(var>U),'or')
    else
        plot(x_outlier(TF),var(TF),'or')
    end
    yline([L U C],":",["Lower Threshold","Upper Threshold","Center Value"])
    if ~isempty(vlineSeparation)
        xline(vlineSeparation,':')
    end
    ylabel(metricLabels{i})
end
xlim([0.5,length(xLabels)+0.5])

xticks(1:length(xLabels))
xticklabels(xLabels)
set(ax,'TickLabelInterpreter', 'none')
linkaxes(ax,'x')

if ~isempty(titlePlot)
    sgtitle(titlePlot)
end

end
