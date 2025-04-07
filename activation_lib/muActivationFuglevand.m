function [muActivation, totalActivation] = muActivationFuglevand(spikeTrains, Ap, tc, initialValues, varargin)
% muActivation calculates motor unit activation profiles.
% totalActivation corresponds to the sum (MU specific activation dynamics)
% Inputs:
% - spikeTrains: Matrix of spike trains (rows: time steps, columns: motor units)
% - Ap: Amplitude parameter
% - tc: Time constant parameter
% - initialValues: Neural activation history (if no history: zeros of size 2 x nMUs)
% Optional inputs (Name-Value pairs):
% - 'fs': Sampling frequency (default: 2048)
% - 'normalization': Type of normalization (default: 0)
% - 'twitchModel': Twitch model type (default: 'Fuglevand')
% - 'twitchSaturation': Enable twitch saturation (default: false)
% - 'satLevel': Saturation level (default: 0.3)
% - 'satSpeed': Saturation speed (default: 100)

%% default values
fs= 2048;
normalization = 0;
% todo improve saturation formulation
numMUs = size(spikeTrains,2);
tetanusSaturation = false;
satLevel = 0.3;
satSpeed = 100;
normSat = numMUs.* satLevel;

%% Process optional input parameters using a loop and switch case
for i = 1:2:length(varargin)
    param = varargin{i};
    value = varargin{i + 1};
    switch param
        case 'fs'
            fs = value;
        case 'normalization'
            normalization = value;
        case 'twitchSaturation'
            twitchSaturation = value;
        case 'satLevel'
            satLevel = value;
        case 'satSpeed'
            satSpeed = value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end

%% initialize parameters
muActivation = [initialValues([size(initialValues,1)-1,size(initialValues,1)],:);
    zeros(size(spikeTrains))];
spikeTrains = [zeros(2,size(spikeTrains,2));spikeTrains];

T = 1/fs;

%% normalization method
if normalization==1 % rec curve witwith generic curve
    updatedTwitchFactor = zeros([1,size(spikeTrains,2)]);
    r=[0,0.4,0.75,1];
    f=[0,0.2,0.5,1]; %normalized filtred force
    recCurve = fit(f',r','poly2');
    normFactor = 1./numMUs; % this goes together with recCurve
elseif normalization ==2  % rec cruve  Duchateau and Enoka 2020 Enokja https://journals.physiology.org/doi/epdf/10.1152/japplphysiol.00290.2021
    % and Oye et al.
    m = 0.045; % from duchateau and enoka 2022, so far I don't know models for other muscles , so going for the TA as generic
    normFactor = 1./numMUs;
    updatedTwitchFactor = zeros([1,size(spikeTrains,2)]);
elseif normalization ==3  %only numMUs
    normFactor = 1./numMUs;
elseif normalization ==4
    normFactor = 1./numMUs;
    Ap = atanh(Ap);
end
%% Start sample per sample decoding
for n = 3:length(spikeTrains)
    
    %if sum(spikeTrains(n,:))>0 % to ensure that only after the first spike compute the models
        if normalization==0   % test og normalization
            muActivation(n,:) = (2*exp(-T./tc))'.*muActivation(n-1,:)  -  (exp((-2*T)./tc))'.*muActivation(n-2,:)  ...
                + (((Ap*T)./tc).*exp(1-T./(tc)))'.*spikeTrains(n-1,:);
        elseif normalization ==1
            % updateTwitchFactor updates the twitch only when a new twitch appears
            updatedTwitchFactor(logical(spikeTrains(n-1,:))) = 1.5*recCurve(abs(force(n)))+0.5;
            muActivation(n,:) = (2*exp(-T./tc))'.*muActivation(n-1,:)  -  (exp((-2*T)./tc))'.*muActivation(n-2,:)  ...
                + ((((normFactor)*Ap*T)./tc).*exp(1-T./(tc)))'.*updatedTwitchFactor.*spikeTrains(n-1,:);
            % practically the same the one above waits for next twitch
            %                     neuralActivation(n,:) = (2*exp(-T./tc))'.*neuralActivation(n-1,:)  -  (exp((-2*T)./tc))'.*neuralActivation(n-2,:)  ...
            %                         + ((((normFactor(i)*recCurve(force(n)))*Ap*T)./tc).*exp(1-T./(tc)))'.*excitation(n-1,:);
        elseif normalization ==2 %normalizing according to recruitmentc coding
            updatedTwitchFactor(logical(spikeTrains(n-1,:))) = ((log(100*abs(force(n-1)))/m) - 2.1178)./100;
            if  updatedTwitchFactor(logical(spikeTrains(n-1,:)))>=1 %saturating to 1
                updatedTwitchFactor(logical(spikeTrains(n-1,:)))=1;
            end
            muActivation(n,:) = (2*exp(-T./tc))'.*muActivation(n-1,:)  -  (exp((-2*T)./tc))'.*muActivation(n-2,:)  ...
                + ((((normFactor)*Ap*T)./tc).*exp(1-T./(tc)))'.*updatedTwitchFactor.*spikeTrains(n-1,:);
        elseif normalization >=3
            muActivation(n,:) = (2*exp(-T./tc))'.*muActivation(n-1,:)  -  (exp((-2*T)./tc))'.*muActivation(n-2,:)  ...
                + (((((normFactor))*abs(force(n)-force(n-1))*Ap*T)./tc).*exp(1-T./(tc)))'.*spikeTrains(n-1,:);
        end
        if tetanusSaturation
            neuralActivationTransf(n,:) = normSat*(1-1*exp(-satSpeed*(neuralActivation(n,:))))./(1+1*exp(-satSpeed*(neuralActivation(n,:))));
        end
    %end
end
if tetanusSaturation
    muActivation = neuralActivationTransf;
end
muActivation(1:2,:)=[]; % elimintes initial values 
totalActivation = sum(muActivation,2);

end
