function [signal,decompParameters] = decompFastICA_v21(EMGfilt,varargin)
%DECOMPFASTICA Perform FastICA decomposition on EMG signals.
%   [SIGNAL, DECOMPPARAMETERS] = DECOMPFASTICA(EMGFILT, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%   performs FastICA decomposition on the input EMG signals (EMGFILT) with optional
%   parameters specified as parameter-value pairs.
%
%   Optional Parameters:
%   - 'fs': Sampling frequency (default: 2048)
%   - 'nbIterations': Number of iterations (default: 50)
%   - 'nbextchan': Number of external channels (default: 1000)
%   - 'SILthreshold': SIL threshold for filtering MU filters (default: 0.85)
%   - 'emgMask': Mask for selecting specific EMG channels (default: all channels)
%   - 'preOptFilters': Pre-optimized filters for initialization (default: [])
%   - 'showPlots': Flag to control whether to show plots or not (default: true)
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
SILthreshold = 0.85;
preOptFilters = [];
showPlots = false;
tanh_denoise = [];
normtanh = [];
flagManualRemoval=0;
removeDuplicates=1;
removeOutliers =1;
cut = 1;
dontremoveanything=0;
optimizeTanh = false;
onlineDecomp = true;

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
        case 'SILthreshold'
            SILthreshold = value;
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
        case 'dontremoveanything'
            dontremoveanything= value;
        case 'normIPT'
            normIPT=value;
        case 'cut'
            cut = value;
        case 'onlineDecomp'
            onlineDecomp = value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end

% Signal Extension
exFactor = round(nbextchan/size(EMGfilt(emgMask,:),1));
eSIG = extend(EMGfilt(emgMask,:),exFactor);


if ~isempty(tanh_denoise)
    normtanh =  tanh_denoise*std(eSIG,[],2);
    eSIG = normtanh.*tanh(eSIG./normtanh);
end

% Signal Whitening
[E, D] = pcaesig(eSIG);
[wSIG, whiteningMatrix, dewhiteningMatrix] = whiteesig(eSIG, E, D);
X = wSIG;
% Removing the edges
if cut > 0
    X = wSIG(:,round(fs):end-round(fs));
end

%X = wSIG(:,1:length(EMGfilt)); % Initialize X (whitened signal), then X: residual
% if ~isempty(tanh_denoise)
%     normtanh =  tanh_denoise*std(X,[],2);
%     if optimizeTanh
%     [normtanh, X] = denoiseTanH(X,1);
%     end
%     X = normtanh.*tanh(X./normtanh);
%     Apply hyperbolic tangent transformation to each row
%     for i = 1:size(X, 1)
%         attenuation_factor = 0.01;
%         Extract the current row
%         data = X(i, :);
% 
%         Apply hyperbolic tangent transformation
%         data_tanh = tanh(data * attenuation_factor);
% 
%         Store the compressed row back into the matrix
%         X(i, :) = data_tanh;
%     end
% end
% %normTmp
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
    SILthreshold = 0.85;
    % TODO make warning or force to include normIPT
end
MUFilters = zeros(size(X,1), nbIterations); % only reliable vectors
SIL = zeros(1,nbIterations);
COV = zeros(1,nbIterations);
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


    maxiter = 100; % max number of iterations for the fixed point algorithm

    
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
        ISI = diff(Distime{j}/fs); % Interspike interval
        COV(j) = std(ISI)/mean(ISI); % Coefficient of variation
        [w, Distime{j},  SIL(j),PulseT(j,:),normIPT(j),centroids(j,:)] = ...
            maximizeSIL(Distime{j}, X, SIL(j),w, fs,showPlots,h,ax);
        %[PulseT(j,:), Distime{j},SIL(j),normIPT(j),centroids(j,:)] = getspikes(MUFilters(:,j), X, fs);
    end
    
    MUFilters(:,j) = w;
    
    tit = ['Ite : ' num2str(j) ' - SIL : ' num2str(round(SIL(j),2)) ' - COV : ' num2str(round(COV(j),2))];
    if ~flagPreOptFilters
        actind(idx1(j)) = 0; % remove the previous vector
    end
    waitbar((j)/(nbIterations),f,tit);
end
% keeps track of good MUs, good for reusing filters
idsnew =  1:nbIterations;

% Filter out MUfilters below the SIL threshold
idsnew(SIL < SILthreshold) = [];
MUFilters(:,SIL < SILthreshold) = [];
COV(SIL < SILthreshold) = [];
PulseT(SIL < SILthreshold,:)=[];
Distime(SIL < SILthreshold)=[];
normIPT(SIL < SILthreshold)=[];
centroids(SIL < SILthreshold,:)=[];
SIL(SIL < SILthreshold) = [];
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

% if onlineDecomp 
%     size(EMGfilt)
%     decompParameters_nonwhit.MUFilters = dewhiteningMatrix*MUFilters;
%     [decompParameters_nonwhit.PulseT,decompParameters_nonwhit.Distime, ~] = getPulseT(EMGfilt, emgMask, decompParameters_nonwhit.MUFilters, fs);
%     % % Get the decomposition parameters for the biofeedback
%     % decompParameters_nonwhit.MUFilters = getMUfilters(EMGfilt, emgMask, Distime);
%     % [decompParameters_nonwhit.PulseT, decompParameters_nonwhit.Distime, ~] = getPulseT(EMGfilt, emgMask, decompParameters_nonwhit.MUFilters, fs);
%     [decompParameters_nonwhit.extensionfactor, decompParameters_nonwhit.iReSIG, decompParameters_nonwhit.norm, decompParameters_nonwhit.centroid] = getonlineparameters(EMGfilt, emgMask, decompParameters_nonwhit.MUFilters, exFactor, fs);
% end 

if showPlots
    % Plotting
    figure;
end
[spikeTrains, meanDR,Distime, ~, ~, ~, ~, removedMUs] = getFiringProperties(Distime,...
        'outlierFlag',removeOutliers,'flagPlotSpikeTrains',showPlots,'flagPlotDR',showPlots,...
        'tVec',linspace(0,length(PulseT)/fs,length(PulseT)),...
        'flagManualRemoval',flagManualRemoval,'dontremoveanything',dontremoveanything);
if ~dontremoveanything
idsnew(removedMUs) = [];
MUFilters(:,removedMUs) = [];
COV(removedMUs) = [];
PulseT(removedMUs,:)=[];
normIPT(removedMUs)=[];
centroids(removedMUs,:)=[];
SIL(removedMUs) = [];
end

signal.spikeTrains = spikeTrains;
signal.PulseT = PulseT;
signal.SIL = SIL;
signal.COV = COV;
signal.Distime= Distime;
signal.idsnew = idsnew;
signal.meanDR =meanDR;
decompParameters.EMGmask = emgMask;
decompParameters.extensionfactor = exFactor;
decompParameters.iReSIG = pinv(eSIG*eSIG'/length(eSIG));
decompParameters.normIPT = normIPT;
decompParameters.normtanh =normtanh;
decompParameters.centroid = centroids;
decompParameters.whitMat = whiteningMatrix;
decompParameters.dewhitMat = dewhiteningMatrix;
decompParameters.MUFilters = MUFilters;

end