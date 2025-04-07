function emg_reRef = re_referencing(emg_data, varargin)
    % INPUT:
    % emg_data -> matrix of EMG signals
    % floating_chnls -> noisy ans silent channels output of the identify_FluctuatingChannels
    % ref_case -> there are three possibile re-referencing possibilities:
    %             1: the reference is just one and is the average of all good channels
    %             2: thre are 8 reference signals, one for each column
    %             3: there are 2 reference signals, one for dorsi-flexors
    %             channels and one for plantar-flexors channels
    % plot_enable -> 0 or 1 if you don't or want to plot a figure
    % 
    %
    % OUTPUT
    % emg_noRef -> re-referencing against the average of all good channels
    % emg_noRef_2 -> re-referencing of each column
    % emg_noRef_3 -> re-referencing using each DF (column 1 and 2) and PF(columns 3 to 8) electrodes ref 
 
% Reshape EMGraw if necessary
flagTransposed = false;
if size(emg_data, 1) < size(emg_data, 2)
    emg_data = emg_data';
    flagTransposed = true;
end

 % Default values for optional parameters
 % usually we have only 2 grids so emgdata is split by 2
referenceConfig =  {[1:round(size(emg_data, 2)/2)],[round(size(emg_data, 2)/2)+1:size(emg_data, 2)]};
emgMask = ones(1,max(max(referenceConfig{:})));
plotFlag = false;

% Process optional input parameters using a loop and switch case
for i = 1:2:length(varargin)
    param = varargin{i};
    value = varargin{i+1};
    switch param
        case 'referenceConfig'
            referenceConfig = value;
        case 'emgMask'
            emgMask = value;
        case 'plotFlag'
            plotFlag = value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end  
emg_reRef = emg_data;
% normF =vecnorm(emg_reRef,2,1);
% emg_reRef = emg_reRef./normF; % works worse by normalizing 

for i=1:length(referenceConfig)
    emgMaskTempSameGrid = emgMask(referenceConfig{i});
    emgGoodChannelsSameGrid = referenceConfig{i}(emgMaskTempSameGrid);
    emgBadChannelsSameGrid = referenceConfig{i}(~emgMaskTempSameGrid);
  
    emg_reRef(:,emgGoodChannelsSameGrid) = emg_reRef(:,emgGoodChannelsSameGrid)-mean(emg_reRef(:,emgGoodChannelsSameGrid),2);
    emg_reRef(:,emgBadChannelsSameGrid) = emg_reRef(:,emgBadChannelsSameGrid)-mean(emg_reRef(:,emgBadChannelsSameGrid),2);
        
    if plotFlag
        sampleChannels= emgGoodChannelsSameGrid(randi(length(emgGoodChannelsSameGrid),[2,1]));
        ax(i,1) = subplot(length(referenceConfig),2,2*(i-1)+1);
        plot((emg_data(:,sampleChannels)))
        ylabel(ax(i,1),['Grid',num2str(i),newline,'EMG'])
        legend(num2str(sampleChannels'))
        ax(i,2) = subplot(length(referenceConfig),2,2*(i-1)+2);
        plot((emg_reRef(:,sampleChannels)))
        legend(num2str(sampleChannels'))
    end
end

% emg_reRef = emg_reRef.*normF; % original magnitude  % works worse by normalizing 

if plotFlag
    % xlabel for the bottom ones (last i)
    xlabel(ax(i,1),'Samples')
    xlabel(ax(i,2),'Samples')
    % title for the top ones
    sgtitle('Re-referencing')
    title(ax(1,1),'Before')
    title(ax(1,2),'After')
    linkaxes(ax,'x')
end
    
    %% Returning to original shape
if flagTransposed
    emg_data = emg_data';
    emg_reRef = emg_reRef';
end

end