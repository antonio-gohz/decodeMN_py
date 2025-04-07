function [signal,decompParameters] = refineMUFilters(EMGs,SILs,decompParameters,varargin)
% simulatedRTdecomp performs real-time decomposition of EMG signals.
% It returns decomposed signals, spike trains, and activation dynamics.
%
% Inputs:
%   EMGs: EMG signals (columns represent muscles, rows represent samples)
%   decompParameters: Decomposition parameters
%   varargin: Optional input parameters
%       - 'fs': Sampling frequency (default: 2048)
%
% Outputs:
%   signal: Decomposed EMG signals and spike trains
%   decompParameters: Decomposition parameters

fs=2048;
showPlots=0;

% Process optional input parameters using a loop and switch case
for i = 1:2:length(varargin)
    param = varargin{i};
    value = varargin{i+1};
    switch param
        case 'fs'
            fs = value;
        case 'showPlots'
            showPlots = value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end

if size(EMGs, 1) > size(EMGs, 2)
    EMGs = EMGs';
end

X = extend(EMGs(decompParameters.EMGmask,:), decompParameters.extensionfactor);
X = X(:,1:size(EMGs, 2));
X = decompParameters.whitMat*X;
X = decompParameters.normtanh.*tanh(X./decompParameters.normtanh);

if showPlots
    %axes properties
    colors =[ 0.2667    0.4667    0.6667;
        0.9333    0.4000    0.4667;0.2667    0.4667    0.6667];
    ylimits = [-0.5,1;-1,1];
    xlimits = [0,size(X, 2)/fs;1,size(decompParameters.MUFilters,1)];
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
signal.spikeTrains= zeros(size(EMGs, 2),size(decompParameters.MUFilters,2));
signal.Pulsetrain= zeros(size(decompParameters.MUFilters,2),size(EMGs, 2));
signal.SIL= zeros(1,size(decompParameters.MUFilters,2));

%f = waitbar(0,'Decomposition...');
for j = 1:size(decompParameters.MUFilters,2)
    %signal.spikeTrains(:,j);
    %ipt = signal.Pulsetrain(j,:);
    SIL = SILs(j);
    w = decompParameters.MUFilters(:,j);
    normIPT = decompParameters.normIPT(j);
    centroids = decompParameters.centroid(j,:);

    [wlast, spikeslast, SILlast,iptlast,normIPTlast,centroidslast] = ...
        maximizeSIL_online(X, SIL,w, fs,showPlots,h,ax,normIPT,centroids);
    
    signal.spikeTrains(:,j) = spikeslast;
    signal.Pulsetrain(j,:) = iptlast;
    signal.SIL(j) = SILlast;

    decompParameters.MUFilters(:,j)=wlast;
    decompParameters.normIPT(j)=normIPTlast;
    decompParameters.centroid(j,:)=centroidslast;
end

