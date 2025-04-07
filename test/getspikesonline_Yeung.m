function  [PulseTout, SpikeTrains, decompParameters, esample2] = getspikesonline_Yeung(EMGtmp,decompParameters, esample2, fs)

% onlineDecomp performs real-time decomposition of EMG signals.
% It returns decomposed signals, spike trains, and activation dynamics.
%dd
% Inputs:
%   EMGs: EMG signals (columns represent muscles, rows represent samples)
%   decompParameters: Decomposition parameters
%   varargin: Optional input parameters
%       - 'fs': Sampling frequency (default: 2048)
%       - 'idxmuscle1': Index of the muscle to analyze (default: 1)
%       - 'refreshRate': Refresh rate (default: 64) [Hz]
%       - 'flagVis': Visualization flag (default: 0)
%       - 'SmoothingEditFieldDR': Smoothing factor for DR (default: 0.2)[s]
%       - 'plotTimeRange': Time range for plotting (default: 10)[s]
%       - 'flagAD': Activation dynamics flag (default: 0)
%       - 'aps': Peak amplitudes for Activation dynamics 
%       - 'tcs': contaction times Activation dynamics [s]
%       - 'stateMask4MUFilters': 
%
% Outputs:
%   signal: Decomposed EMG signals and spike trains
%   decompParameters: Decomposition parameters
%   activation: Activation dynamics

% Check if EMGs needs to be transposed
lambda = 0.985;
zvalue = 3.3;
MUFilters = decompParameters.MUFilters;
%eSIG = extend(EMGtmp,decompParameters.extensionfactor);
%Esample 1 for which the train has to be extracted.
%esample1 = eSIG(:, 1:size(EMGtmp, 2));
esample1 = EMGtmp;
% fill zeros at the beginning of the extended data with data from the previous block (esample2)
%esample1(:, 1:decompParameters.extensionfactor-1) = esample1(:, 1:decompParameters.extensionfactor-1) + esample2;
% retrieve the extend2 data for the next block
%esample2 = eSIG(:, size(EMGtmp, 2) + 1:end);


ReSIG_win = esample1*esample1'/length(esample1);
ReSIG = decompParameters.ReSIG*lambda + (1-lambda)*ReSIG_win; %(1-lambda)
iReSIGt = pinv(ReSIG);
buffer_size = size(decompParameters.buffer_eSIG,2);

for i = 1:size(MUFilters,2)
    PulseTtmp = (MUFilters(:,i)'*iReSIGt)*esample1; % '*iReSIGt)
    % PulseTtmp = PulseTtmp(1:size(EMG,2));
    PulseT = PulseTtmp .* abs(PulseTtmp);
    [~,spikes] = findpeaks(PulseT, 'MinPeakDistance', round(fs*0.02)); % 4b: Peak detection
    spikesTmp = spikes;
    %spikesTmp(isoutlier(PulseT(spikes),'percentiles',[0,99]))=[]; % NEW add remove the outliers
    PulseT= tanh(PulseT/decompParameters.norm(i));
        PulseZ = zscore(PulseT);
    % See if z-score is above treshold.
    PotSpikes = find(PulseZ(spikesTmp) > zvalue);
    peaksnoise = setdiff(spikesTmp, spikesTmp(PotSpikes));
    PotSpikes = PotSpikes(PulseT(spikesTmp(PotSpikes)) > 0.10);
    if ~isempty(PotSpikes) 
        buffer_eSIG = cat(2,decompParameters.buffer_eSIG(:,:,i),esample1(:,spikesTmp(PotSpikes)));
        buffer_spikes = cat(2,decompParameters.buffer_spikes(i,:), PulseT(spikes(PotSpikes)));
        decompParameters.buffer_eSIG(:,:,i) = buffer_eSIG(:,end-buffer_size+1:end);
        decompParameters.buffer_spikes(i,:) = buffer_spikes(end-buffer_size+1:end);
        % For now the centroids, same updates as are proposed in the Yeung
        % et al. algorithm. However, this time, we use the average
        % amplitude of the normalized pulse trains. 
        decompParameters.centroid(i,1) = mean(decompParameters.buffer_spikes(i,:), 'omitnan');
        decompParameters.ReSIG = ReSIG;
        if ~isempty(peaksnoise)
            decompParameters.centroid(i,2) = lambda*decompParameters.centroid(i,2)+ (1-lambda)*mean(PulseT(peaksnoise));
        end
        decompParameters.MUFilters(:,i) = mean(decompParameters.buffer_eSIG(:,:,i),2, 'omitnan');
    end
    %Classification of the peaks after updates of the centroids:
    PulseTout(i,:) = PulseT;
    %[spikes1, ~] = islocalmax(PulseT, 1, 'MinSeparation', round(fs*0.02));
    SpikeTrains(i,:) = zeros(size(PulseT'));
    SpikeTrains(i,spikesTmp(PotSpikes)) = (abs(PulseT(spikesTmp(PotSpikes))'  - decompParameters.centroid(i,2)) > abs(PulseT(spikesTmp(PotSpikes))' - decompParameters.centroid(i,1)));
end


