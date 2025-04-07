function [signal,decompParameters] = decompFastICA_v2(EMGfilt,varargin)
%DECOMPFASTICA Perform FastICA decomposition on EMG signals.
%   [SIGNAL, DECOMPPARAMETERS] = DECOMPFASTICA(EMGFILT, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%   performs FastICA decomposition on the input EMG signals (EMGFILT) with optional
%   parameters specified as parameter-value pairs.
%
%   Optional Parameters:
%   - 'fs': Sampling frequency (default: 2048)
%   - 'nbIterations': Number of iterations (default: 50)
%   - 'nbextchan': Number of extended channels (default: 1000)
%   - 'SILthreshold': SIL threshold for filtering MU filters (default: 0.85)
%   - 'emgMask': Mask for selecting specific EMG channels (default: all channels)
%   - 'preOptFilters': Pre-optimized filters for initialization (default: [])
%   - 'showPlots': Flag to control whether to show plots or not (default: true)
%   - 'refineStrategy': Refinement strategy, COV is recommended for constatn force contractions (default: SIL)
%   - 'peeloff_flag': Peel off flag,  (default: false)
%
%   Output:
%   - SIGNAL: Struct containing discharge times, pulse train, SIL, COV, and more.
%   - DECOMPPARAMETERS: Struct containing decomposition parameters for biofeedback.
%   EXAMPLES:
%

% Default values for optional parameters
fs = 2048;
nbIterations = 50;
nbextchan = 1000;
maxiter = 100; % max number of iterations for the fixed point algorithm
buffersize = 100;
preOptFilters = [];
showPlots = false;
tanh_denoise = [];
normtanh = [];
flagManualRemoval=0;
removeDuplicates=1;
removeOutliers =1;
keepAllMUs=0;
optimizeTanh = false;
peeloff_flag= false;
qc_threshold = 0;
refineStrategy = 'SIL';

if size(EMGfilt, 1) > size(EMGfilt, 2)
    EMGfilt = EMGfilt';
end
emgMask = true(1,size(EMGfilt,1));
% Process optional input parameters using a loop and switch case
for i = 1:2:length(varargin)
    param = varargin{i};
    value = varargin{i+1};
    switch param
        case 'fs'
            fs = value;
        case 'nbIterations'
            nbIterations = value;
        case 'nbextchan'
            nbextchan = value;
        case 'maxiter'
            maxiter = value;
        case 'buffersize'
            buffersize = value;
        case 'qc_threshold'
            qc_threshold = value;
        case 'emgMask'
            emgMask = value;
        case 'preOptFilters'
            preOptFilters = value;
        case 'showPlots'
            showPlots = value;
        case 'tanh_denoise'
            tanh_denoise = value; % number of sigma for sat
        case 'optimizeTanh' % todo add initial value
            optimizeTanh  = value; 
        case 'flagManualRemoval'
            flagManualRemoval =value;
        case 'removeOutliers'
            removeOutliers =value;
        case 'removeDuplicates'
            removeDuplicates =value;
        case 'keepAllMUs' % previously 'dontremoveanything'
            keepAllMUs= value;
        case 'normIPT'
            normIPT=value;
        case 'refineStrategy'
            refineStrategy=value;
        case 'peeloff_flag'
            peeloff_flag=value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end

if qc_threshold==0
    if strcmp(refineStrategy,'SIL')
        qc_threshold = 0.85;
    elseif strcmp(refineStrategy,'COV')
        qc_threshold = 0.6;
    end
end

% Signal Extension
exFactor = round(nbextchan/size(EMGfilt(emgMask,:),1));
eSIG = extend(EMGfilt(emgMask,:),exFactor);

% Signal Whitening
[E, D] = pcaesig(eSIG);
[wSIG, whiteningMatrix, dewhiteningMatrix] = whiteesig(eSIG, E, D);

% Removing the edges
%X = wSIG(:,round(fs):end-round(fs));

X = wSIG(:,1:length(EMGfilt)); % Initialize X (whitened signal), then X: residual
if ~isempty(tanh_denoise)
    normtanh =  tanh_denoise*std(X,[],2);
    if optimizeTanh
    [normtanh, X] = denoiseTanH(X,1);
    end
X = normtanh.*tanh(X./normtanh);
end
%normTmp
% normTmp = maxk(X,10,2);
% normTmp =  mean(normTmp,2);
% baselineTmp = maxk(X(:,1650:1800),10,2);
% baselineTmp =  mean(baselineTmp,2);
% if mean(normTmp-baselineTmp)<2
%     baselineTmp =  zeros(size(baselineTmp));
% end
% %signsOG = sign(X);
% X = tanh(2*X./(normTmp-baselineTmp));
%X = X.*normTmp;
% FastICA
% Initialize matrix B (n x m) n: separation vectors, m: iterations
% Initialize matrix MUFilters to only save the reliable filters
% Intialize SIL and PNR
% check pre-opt filters and update nbIterations
if isempty(preOptFilters)
    flagPreOptFilters =false;
    Xtmp = X;
    Xtmp(isoutlier(X)) = 0; % remove artifacts from activity index
    actind = sum(abs(Xtmp).^2,1);    
    % Find the index where the square of the summed whitened vectors is
    % maximized and initialize W with the whitened observations at this time
    idx1 = zeros(1,nbIterations);
    normIPT= zeros(1,nbIterations);
else
    flagPreOptFilters =true;
    nbIterations = size(preOptFilters,2);
    qc_threshold = 0.85;
    % TODO make warning or force to include normIPT
end
MUFilters = zeros(size(X,1), nbIterations); % only reliable vectors
SIL = zeros(1,nbIterations);
COV = zeros(1,nbIterations);
qc_metric = zeros(1,nbIterations);%quality control metric
PulseT= zeros( nbIterations,size(X,2));
Distime= cell(1,nbIterations); 
centroids = zeros( nbIterations,2);

if showPlots
    %axes properties
    colors =[ 0.2667    0.4667    0.6667;
        0.9333    0.4000    0.4667;0.2667    0.4667    0.6667];
    ylimits = [-0.5,1;-1,1];
    xlimits = [0,length(X)/fs;1,length(MUFilters)];
    ylabels = {'IPTs','MU filter'};
    lineStyles = {'-','none','-'};
    markerStyles = {'none','o','none'};
    axesLineConfig = [1,1;    % subplot 1 -signal 1 (IPT)
        1,2;       % subplot 1 -signal 2 (IPT peaks)
        2,3];      % % subplot 2 -signal 3 (MU filter)
    [h,ax] = livePlots.preparePlots(1,'linewidth',1,'colors',colors,...
        'ylimits',ylimits ,'xlimits',xlimits,'ylabels',ylabels,...
        'lineStyles',lineStyles,'markerStyles',markerStyles,...
        'axesLineConfig',axesLineConfig);
else
    h=[]; ax=[];
end



    
f = waitbar(0,'Decomposition...');
for j = 1:nbIterations
    
    if ~flagPreOptFilters
        [~, idx1(j)] = max(actind);
        w = X(:, idx1(j)); % Initialize w
        w = w - MUFilters * MUFilters' * w; % Orthogonalization
        w = w / norm(w); % Normalization  
        %       3a: Fixed point algorithm (end when sparsness is maximized)
        w = fixedpointalg(w, X, MUFilters , maxiter, 'logcosh',fs,showPlots,h,ax);
    else
        w = preOptFilters(:,j); % Initialize w
        %w = sum(X(:,preOptFilters{j}),2); % update W by summing the spikes
    end    
    
    % Step 4: Maximization of SIL:  Initialize SIL Step 4a => 4e
    [PulseT(j,:), Distime{j},SIL(j),normIPT(j),centroids(j,:)] = getspikes(w, X, fs);
    
    if length(Distime{j}) > 10 % If the number of peaks is low, skip the second algorithm
        if strcmp(refineStrategy,'SIL')
            [w, Distime{j},  SIL(j),PulseT(j,:),normIPT(j),centroids(j,:)] = ...
                maximizeSIL(Distime{j}, X, SIL(j),w, fs,showPlots,h,ax);
            ISI = diff(Distime{j}/fs); % Interspike interval
            COV(j) = std(ISI)/mean(ISI); % Coefficient of variation

            qc_metric(j) = SIL(j);

        elseif strcmp(refineStrategy,'COV')
            ISI = diff(Distime{j}/fs); % Interspike interval
            COV(j) = std(ISI)/mean(ISI); % Coefficient of variation
            [w,  Distime{j}, COV(j),SIL(j),PulseT(j,:),normIPT(j),centroids(j,:)] = ...
                minimizeCOVISI(w, X, COV(j), fs,showPlots,h,ax);
            [~, ~, SIL(j)] = calcSIL(X, w, fs);

            qc_metric(j) = COV(j);
        end

        if peeloff_flag % Peel-off of the (reliable) source
            X = peeloff(X, spikes, round(str2double(app.FrequencyDropDown.Value)), 0.025);
        end

    end
    % whitened MU filters
    MUFilters(:,j) = w;
    
    
    tit = ['Ite : ' num2str(j) ' - SIL : ' num2str(round(SIL(j),2)) ' - COV : ' num2str(round(COV(j),2))];
    if ~flagPreOptFilters
        actind(idx1(j)) = 0; % remove the previous vector
    end
    waitbar((j)/(nbIterations),f,tit);
end
% keeps track of good MUs, good for reusing filters
idsnew =  1:nbIterations;

% Filter out MUfilters below the SIL threshold or above COV threshold
%changes sign accordingly (i.e. <, > )
if strcmp(refineStrategy,'SIL')
    qc_sign= 1;
elseif strcmp(refineStrategy,'COV')
    qc_sign = -1;
end
idsnew(qc_sign*qc_metric < qc_sign*qc_threshold) = [];
MUFilters(:,qc_sign*qc_metric < qc_sign*qc_threshold) = [];
COV(qc_sign*qc_metric < qc_sign*qc_threshold) = [];
PulseT(qc_sign*qc_metric < qc_sign*qc_threshold,:)=[];
Distime(qc_sign*qc_metric < qc_sign*qc_threshold)=[];
normIPT(qc_sign*qc_metric < qc_sign*qc_threshold)=[];
centroids(qc_sign*qc_metric < qc_sign*qc_threshold,:)=[];
SIL(qc_sign*qc_metric < qc_sign*qc_threshold) = [];
%MUFilters(:,COV > COVthreshold) = [];
%MUFilters = dewhiteningMatrix * MUFilters;

% Get the pulse train for the entire signal
%[PulseT, Distime, ~] = getPulseT(EMGfilt, ~emgMask, MUFilters, fs,nbextchan);

% Remove duplicates
close(f)
if removeDuplicates
    [PulseT, Distime,idsnewVec] = remduplicates(PulseT, Distime, Distime, round(fs/40), 0.00025, 0.3, fs);
    SIL = SIL(idsnewVec);
    COV = COV(idsnewVec);
    MUFilters = MUFilters(:,idsnewVec);
    normIPT=normIPT(idsnewVec);
    centroids=centroids(idsnewVec,:);
    idsnew = idsnew(idsnewVec);
end

if showPlots
    % Plotting
    figure;
end
[spikeTrains, dischargeRates,Distime, ~, ~, ~, ~, removedMUs] = getFiringProperties(Distime,...
        'outlierFlag',removeOutliers,'flagPlotSpikeTrains',showPlots,'flagPlotDR',showPlots,...
        'tVec',linspace(0,length(PulseT)/fs,length(PulseT)),'fs_MN',fs,...
        'flagManualRemoval',flagManualRemoval,'keepAllMUs',keepAllMUs);
if ~isempty(removedMUs)
idsnew(removedMUs) = [];
MUFilters(:,removedMUs) = [];
COV(removedMUs) = [];
PulseT(removedMUs,:)=[];
normIPT(removedMUs)=[];
centroids(removedMUs,:)=[];
SIL(removedMUs) = [];
end

% Non-whitened MU filters
MUFilters = getMUfilters(eSIG, Distime);
% [PulseT, Distime, ~] = getPulseT(EMGfilt, ~emgMask, decompParameters.MUFilters, fs, nbextchan);
[ReSIG, iReSIGt, normIPT, centroids] = getonlineparameters(eSIG, MUFilters, fs);

signal.spikeTrains = spikeTrains';
signal.Pulsetrain = PulseT;
signal.SIL = SIL;
signal.COV = COV;
signal.Dischargetimes= Distime;
signal.idsnew = idsnew;
signal.dischargeRates =dischargeRates';
decompParameters.buffersize = buffersize;  % todo make it an input optional
decompParameters.buffer_eSIG = NaN(size(eSIG,1), buffersize, size(MUFilters,2));
decompParameters.buffer_spikes = NaN(size(MUFilters,2), buffersize);
decompParameters.EMGmask = emgMask;
decompParameters.extensionfactor = exFactor;
decompParameters.ReSIG = ReSIG;
decompParameters.iReSIGt = iReSIGt;
decompParameters.normIPT = normIPT;
decompParameters.normtanh =normtanh;
decompParameters.centroid = centroids;
decompParameters.whitMat = whiteningMatrix;
decompParameters.dewhitMat = dewhiteningMatrix;
decompParameters.MUFilters = MUFilters;

end


