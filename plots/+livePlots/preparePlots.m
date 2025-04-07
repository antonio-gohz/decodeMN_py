function [h,ax,fig] = preparePlots(numMUs,varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%axes properties
ylimits = [0,40;0.5 numMUs+0.5;0,1];
ylabels = {'Discharge rate (Hz)','MU Index','MU-activations (a.u.)'}; 
xlabels = repelem({'Samples'},length(ylabels));
xlimits = [0,10];
lineStyles = {'none','-','-','-'};
markerStyles = {'o','none','none','none'};
colors =[

    0.2667    0.4667    0.6667
    0.9333    0.4000    0.4667
    0.1333    0.5333    0.2000
    0.8000    0.7333    0.2667
    0.4000    0.8000    0.9333
    0.6667    0.2000    0.4667
    0.7333    0.7333    0.7333];
% [axes, line]; DR usually has too lines 'o' for the tip and '-' for the rest'
% default does not include 3rd axes e.g. for Activation dynamics plot [3,4]
axesLineConfig = [1,1;
    1,2;
    2,3];
linewidth=2;
for i = 1:2:length(varargin)
    param = varargin{i};
    value = varargin{i + 1};
    switch param
        case 'axesLineConfig'
            axesLineConfig = value;
        case 'lineStyles'
            lineStyles = value;
        case 'markerStyles'
            markerStyles = value;
        case 'ylimits'
            ylimits = value;
        case 'xlimits'
            xlimits = value;
        case 'ylabels'
            ylabels = value;
        case 'xlabels'
            xlabels = value;
        case 'colors'
            colors = value;
        case 'linewidth'
            linewidth = value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end

fig = figure('units','normalized','outerposition',[0 0 0.7 1]);
%VisualizationBiofeedback.Position = [100 100  690 1000];
fig.Name = 'MU Feddback';

% number of axes to be plot 
for i=1:max(axesLineConfig(:,1)) 
    ax(i) = subplot(max(axesLineConfig(:,1)),1,i);
    ylabel(ax(i), ylabels{i})
    xlabel(ax(i), xlabels{i})    
    ax(i).YLim = ylimits(i,:); %max(decompParameters.DR{idxmuscle1})*1.2];
    ax(i).XLim = xlimits(i,:);
end


% Axes properties
%axes_spikes = axes(VisualizationBiofeedback,'units','normalized','Position',[0.1300    0.4096    0.7750    0.2157]);
%xlabel(ax, 'Time (s)')
% set(ax,{'Color'},{'none'})
set(ax,{'GridColor'},{'none'})
set(ax,{'MinorGridColor'},{[0 0 0]})
set(ax,{'FontSize'},{14})
set(ax,{'Interruptible'},{'off'})
%set(ax,{'HitTest'},{'off'})
%set(ax,{'PickableParts'},{'none'})
set(ax,{'Units'},{'normalized'})

%set(ax,{'Xlim'},{xlimits})
%axes_spikes.Toolbar.Visible = 'off';
%axes_spikes.Position = [0.1 0.1 0.35 0.85];
% number of line sets (h) to be plot 

for i=1:max(axesLineConfig(:,2)) 
    h(:,i)=line(ax(axesLineConfig(i,1)),nan,repelem(nan,numMUs),'Color',colors(i,:),'LineStyle',lineStyles{i},'Marker',markerStyles{i}, 'LineWidth', linewidth, 'MarkerSize', 20);
end

end