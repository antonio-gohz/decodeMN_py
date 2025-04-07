function [signal,decompParametersOut] = pseudoAdaptDecomp(EMGs,decompParameters,varargin)
% simulatedRTdecomp performs real-time decomposition of EMG signals.
% It returns decomposed signals, spike trains, and activation dynamics.
%
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
if size(EMGs, 1) > size(EMGs, 2)
    EMGs = EMGs';
end
size(EMGs)
% default parameters
fs=2048;
idxmuscle1=1;
refreshRate = 10 ; % 64 (15.6 ms) - 10 (100 ms)

flagAD = 0;
aps = [];
tcs = [];
tanh_denoise =[];

flagVis =0;
SmoothingEditFieldDR=0.5; % seconds
plotTimeRange = 10; % seconds
axesLineConfig = [1,1;
    1,2;
    2,3];
stateMask4MUFilters = [];

for i = 1:2:length(varargin)
    param = varargin{i};
    value = varargin{i + 1};
    switch param
        case 'fs'
            fs = value;
        %case 'idxmuscle1'
        %    idxmuscle1 = value;
        case 'refreshRate'
            refreshRate = value;
        %case 'flagAD'
        %    flagAD = value;
        case 'flagVis'
            flagVis = value;
        case 'SmoothingEditFieldDR'
            SmoothingEditFieldDR = value;
        case 'plotTimeRange'
            plotTimeRange = value;
        case 'aps'
            aps = value;
        case 'tanh_denoise'
            tanh_denoise = value;
        case 'tcs'
            tcs = value;
        case 'stateMask4MUFilters'
            stateMask4MUFilters = value;
        case 'tanhNormEMG'
            tanhNormEMG = value;
        case 'tanhNormWEMG'
            tanhNormWEMG = value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end


% Simulating adquisition parameters
eSIG = extend(EMGs(decompParameters.EMGmask,:), decompParameters.extensionfactor);
% truncate
eSIG = eSIG(:,round(fs):end-round(fs));

if ~isempty(tanh_denoise)
    normtanh =  tanh_denoise*std(eSIG,[],2);
    eSIG = normtanh.*tanh(eSIG./normtanh);
end

ComParameters.nsamp = floor(fs/refreshRate); %2048 /8
OnlineParameters.durationonline = length(eSIG);
nwin = floor(OnlineParameters.durationonline/ComParameters.nsamp); % number of windows Duration/ fs / refrates(refreshrate)

% Pre allocate empty matrices for decomposition
EMG = zeros(size(EMGs));
time_DR = linspace(0,size(EMG,2)/fs, nwin);
time_spikes_AD = linspace(0,size(EMG,2)/fs, size(EMG,2));
EMGtmp = zeros(size(EMG,1),ComParameters.nsamp);
esample2 = zeros((size(EMG,1)-sum(~decompParameters.EMGmask))*decompParameters.extensionfactor(idxmuscle1),decompParameters.extensionfactor(idxmuscle1)-1);
spikeTrains = zeros(nwin*ComParameters.nsamp, size(decompParameters.MUFilters,2));
signal.Pulsetrain = zeros(size(decompParameters.MUFilters,2), ComParameters.nsamp*nwin);
noise_centroids = ones(ComParameters.nsamp,size(decompParameters.MUFilters,2)) .* (decompParameters.centroid(:,1)');
spike_centroids = ones(ComParameters.nsamp,size(decompParameters.MUFilters,2)) .* (decompParameters.centroid(:,2)');



DR = zeros(nwin, size(decompParameters.MUFilters,2));
ct_dec = zeros([1,nwin]);


if flagVis
    SmoothingEditFieldDR = round(SmoothingEditFieldDR*refreshRate);
    plotLengthWindows =  round(refreshRate*plotTimeRange/2);
    numMUs = size(decompParameters.MUFilters,2);
    [h,ax,~]=preparePlots(numMUs,'axesLineConfig',axesLineConfig);
end

k = 1;
f = waitbar(0,'Decomposition...');

while k <= nwin
    EMG1(:,(k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp) = eSIG(:,(k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp );
    EMGtmp = EMG1(:, (k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp);
    % EMG(:,(k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp) = EMGs(:,(k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp );
    % EMGtmp = EMG(decompParameters.EMGmask, (k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp);

    %EMGtmp = tanhNormEMG.*tanh(EMG(decompParameters.EMGmask, (k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp)./tanhNormEMG);
    tic
%     [signal.Pulsetrain(:,(k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp), spikeTrains((k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp, :), esample2] = ...
%         getspikesonline(EMGtmp, decompParameters.extensionfactor(idxmuscle1), esample2, MUfilt, decompParameters.norm, noise_centroids, spike_centroids, ComParameters.nsamp, fs);
%     [signal.Pulsetrain(:,(k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp), spikeTrains((k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp, :), esample2] = ...
%         getspikesonline(EMGtmp, decompParameters.extensionfactor(idxmuscle1), esample2, decompParameters.MUFilters, decompParameters.iReSIG, decompParameters.norm, noise_centroids, spike_centroids, ComParameters.nsamp, fs,tanhNormWEMG);
%[signal.Pulsetrain(:,(k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp), spikeTrains((k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp, :), esample2] = ...
%       getspikesonline(EMGtmp, decompParameters.extensionfactor(idxmuscle1), esample2, decompParameters.MUFilters,  decompParameters.whitMat, decompParameters.normIPT, decompParameters.centroid,noise_centroids, spike_centroids, ComParameters.nsamp, fs,decompParameters.normtanh);
    [signal.PulseT(:,(k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp), spikeTrains((k-1)*ComParameters.nsamp+1:k*ComParameters.nsamp, :), decompParameters, esample2] = ...
        getspikesonline_Yeung(EMGtmp, decompParameters, esample2,fs);
    ct_dec(k)=toc;
     
    % if flagVis
    %     if k>SmoothingEditFieldDR % compute mean of window (non zeroes only)
    %         DR(k,:) = sum(spikeTrains((k-SmoothingEditFieldDR)*ComParameters.nsamp:k*ComParameters.nsamp, :))*fs/(ComParameters.nsamp*SmoothingEditFieldDR);
    %         DR(k,:) = mean(DR(k-round(SmoothingEditFieldDR/3):k,:), 1);
    %     end
    %     if k < plotLengthWindows % if less than half window length
    %         updatePlot({time_DR(k),DR(k,:);time_DR(1:k),DR(1:k,:);
    %             time_spikes_AD(1:k*ComParameters.nsamp),...
    %             (0.8*(spikeTrains(1:k*ComParameters.nsamp,:))+(1:numMUs)-0.4);
    %             time_spikes_AD(1:k*ComParameters.nsamp),...
    %             muActivation(1:k*ComParameters.nsamp,:)/7.5},...
    %             h,ax,[time_DR(1) time_DR(plotLengthWindows*2)],[0,1,1]);
    %     else
    %         updatePlot({time_DR(k),DR(k,:); time_DR(k-plotLengthWindows+1:k),DR(k-plotLengthWindows+1:k,:);
    %             time_spikes_AD((k-plotLengthWindows+1)*ComParameters.nsamp+1:k*ComParameters.nsamp),...
    %             (0.8*(spikeTrains((k-plotLengthWindows+1)*ComParameters.nsamp+1:k*ComParameters.nsamp,:))+(1:numMUs)-0.4);
    %             time_spikes_AD((k-plotLengthWindows+1)*ComParameters.nsamp+1:k*ComParameters.nsamp),...
    %             muActivation((k-plotLengthWindows+1)*ComParameters.nsamp+1:k*ComParameters.nsamp,:)/7.5},...
    %             h,ax,[time_DR(k-plotLengthWindows+1) time_DR(k)+plotTimeRange/2],[0,1,1]); %[h,ax] =
    %     end
    %     pause(ComParameters.nsamp/fs)
    % end
    tit = ['Block : ' num2str(k) , ' of ', num2str(nwin)];
    waitbar((k)/(nwin),f,tit);

    k = k+1;
end

for m = 1:size(spikeTrains,2)
    signal.Distime{m} = find(spikeTrains(:,m)' > 0);
end

% Store signal and activation variables 
signal.spikeTrains = spikeTrains;
signal.ct = ct_dec;
decompParametersOut = decompParameters;
% activation.muActivation= muActivation ;
% activation.totalActivation = totalActivation;
% activation.ct= ct_AD;
close(f)
end