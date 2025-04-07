function [OnlineDecompParameters] = getonlineparameters(eSIG, parameters, fsamp)
%   getonlineparameters: extracts the centroids, norm and buffers from
%   offline decomposition data.
%
%   Input:
%   'eSIG': extended whitened observations, same size as in offline
%   decomposition
%   'parameters': Contains MUfilters, extensionfactor, buffersize
%   'fsamp': Sample frequency
%
%   Output:
%   'OnlineDecompParameters': Contains buffers, MUfilters, centroids, norms, ReSIG and
%   extensionfactor
   


MUFilters = parameters.MUFilters;
%Retrieve parameters from parameters
%
OnlineDecompParameters.iReSIGt = parameters.iReSIGt;
OnlineDecompParameters.ReSIG= parameters.ReSIG;
OnlineDecompParameters.extensionfactor = parameters.extensionfactor;
OnlineDecompParameters.buffer_eSIG = NaN(size(eSIG,1), parameters.buffersize, size(MUFilters,2));
OnlineDecompParameters.buffer_spikes = NaN(size(MUFilters,2), parameters.buffersize);


for i = 1:size(MUFilters,2)
    %CKC algorithm
    PulseTtmp = (MUFilters(:,i)'*OnlineDecompParameters.iReSIGt)*eSIG; 
    % PulseTtmp = PulseTtmp(1:size(EMG,2));
    % calculate peaks (keep negatives)
    PulseT = PulseTtmp .* abs(PulseTtmp);
    %Detect peak with max 50 Hz discharge rate
    [~,spikes] = findpeaks(PulseT, 'MinPeakDistance', round(fsamp*0.02)); % 4b: Peak detection
    spikesTmp = spikes;

    % Exclude really high peaks to avoid bias of MUfilters
    spikesTmp(isoutlier(PulseT(spikes),'percentiles',[0,99]))=[]; % NEW add remove the outliers

    % Calculate norm as mean of 10 highest peaks (MUST BE CHANGED FOR ADAPTIVE ALGORITHM)
    OnlineDecompParameters.norm(i) = mean(maxk(PulseT(spikesTmp),10));
    
    % Normalize PulseTrains and use tanh to avoid high peaks
    PulseT= tanh(PulseT/OnlineDecompParameters.norm(i));

    if length(spikes)>2
        [L,C] = kmeans(PulseT(spikes)', 2); % 4c: Kmean classification
        [~, idx] = max(C); % Spikes should be in the class with the highest centroid
        peakspikes = spikes(L==idx);
        peaksnoise = setdiff(spikes, peakspikes);
        %peakspikes(PulseT(peakspikes)>mean(PulseT(peakspikes))+3*std(PulseT(peakspikes))) = []; % remove the outliers of the pulse train
           

        if parameters.buffersize < length(peakspikes)
            % if the number of peaks is larger than the buffer, include only
            % the last #buffersize peaks in the buffer.

            % Buffer for calculation of the MU filter 
            OnlineDecompParameters.buffer_eSIG(:,:,i) = eSIG(:,peakspikes(end-parameters.buffersize+1:end));

            % Buffer for calculation of the spike centroid (mean of
            % Pulsetrain values on detected peaks)
            OnlineDecompParameters.buffer_spikes(i,:) = PulseT(peakspikes(end-parameters.buffersize+1:end));
            % Spike centroid calulation 
            OnlineDecompParameters.centroid(i,1) = mean(PulseT(peakspikes(end-parameters.buffersize+1:end)));
            % Noise centroid calulation 
            OnlineDecompParameters.centroid(i,2) = mean(PulseT(peaksnoise));
            % Recalculation of MUfilters
            OnlineDecompParameters.MUFilters(:,i) = mean(OnlineDecompParameters.buffer_eSIG(:,:,i), 2);
        else
            % if the number of peaks is smaller than the buffer, include
            % all peaks, but keep buffer the same size 
            
            % Buffer for calculation of the MU filter 
            OnlineDecompParameters.buffer_eSIG(:,end-length(peakspikes)+1:end,i) = eSIG(:,peakspikes);

            % Buffer for calculation of the spike centroid (mean of
            % Pulsetrain values on discharge times)
            OnlineDecompParameters.buffer_spikes(i,end-length(peakspikes)+1:end) = PulseT(peakspikes);
            % Spike centroid calulation 
            OnlineDecompParameters.centroid(i,1) = mean(PulseT(peakspikes),'omitnan');
            % Noise centroid calulation 
            OnlineDecompParameters.centroid(i,2) = mean(PulseT(peaksnoise),'omitnan');

            %Only take into account non-NAN values for the recalculation of
            %the MUfilters as buffer is not full. 
            OnlineDecompParameters.MUFilters(:,i) = mean(OnlineDecompParameters.buffer_eSIG(:,:,i), 2, 'omitnan');
        end
    end

end
