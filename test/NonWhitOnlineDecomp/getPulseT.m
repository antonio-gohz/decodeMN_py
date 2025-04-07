function [PulseT, Distime, exFactor] = getPulseT(eSIG, parameters, fsamp)
%   getPulseT: calculates the pulse trains and dischargetimes based on
%   offline parameters
%
%   Input:
%   'eSIG': extended whitened observations, same size as in offline
%   decomposition
%   'parameters': Contains MUfilters, extensionfactor, buffersize
%   'fsamp': Sample frequency
%
%   Output:
%   'PulseT': matrix with all the pulse trains of each motor unit in each
%   row
%   'Distime': cell with dischargetimes of each motor unit
%   'exFactor': Extension factor for the observations


MUFilters = parameters.MUFilters;
iReSIGt = parameters.iReSIGt;
exFactor = parameters.extensionfactor;

eSIG(1,1:10)
for i = 1:size(MUFilters,2)
    %CKC algorithm
    PulseTtmp = (MUFilters(:,i)'*iReSIGt)*eSIG; % '*iReSIGt)
    % PulseTtmp = PulseTtmp(1:size(EMG,2));
    % calculate peaks (keep negatives)
    PulseT(i,:) = PulseTtmp .* abs(PulseTtmp);
    %find peaks with maximum discharge rate of 50 Hz.
    [~,spikes] = findpeaks(PulseT(i,:), 'MinPeakDistance', round(fsamp*0.02)); % 4b: Peak detection
    spikesTmp = spikes;

    % Exclude really high peaks
    spikesTmp(isoutlier(PulseT(i,spikes),'percentiles',[0,99]))=[]; % NEW add remove the outliers

    % Normalize pulse trains
    normIPT = mean(maxk(PulseT(i,spikesTmp),10));
    PulseT(i,:)= tanh(PulseT(i,:)/normIPT);
    if length(spikes)>2
        [L,C] = kmeans(PulseT(i,spikes)',2); % 4c: Kmean classification
        [~, idx] = max(C); % Spikes should be in the class with the highest centroid
        
        % Distimes are peaks that are in the class of the highest centroid
        Distime{i} = spikes(L==idx);
        % exclude outliers as discharge times.
        Distime{i}(PulseT(i,Distime{i})>mean(PulseT(i,Distime{i}))+3*std(PulseT(i,Distime{i}))) = []; % remove the outliers of the pulse train
    end
end
