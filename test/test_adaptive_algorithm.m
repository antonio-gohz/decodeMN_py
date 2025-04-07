%% LOAD DATA FROM FILES 
% These files contains outputs of the decompFastICA_v21 called
%'signal{type}' with type either train or test and decompParameters
% They also contain OnlineDecompParameters_{type}, which contain the non whitened
% versions of MUFilters, norms, buffers and centroids of the offline
% decomposition that can be used for online decomposition in the adaptive
% algorithm.
% adapt_data.mat contains AdaptDecompParameters_nonwhit and signalAdapt in
% which the trained MU filters are applied to the test data with the
% adaptive algorithm.
clear all
addpath(genpath('./'))

rootData = ['../../Decomposition project - Nathan/Data/BioRob (preliminary study)/Subject0_test Michele/'];
fs = 2048;
load([rootData, 'train_test_data/adapt_data.mat']);
load([rootData, 'train_test_data/train_data.mat']);
load([rootData, 'train_test_data/test_data.mat']);

%% Load the Saved MU filter data
clear all
close all
% overkill adds all libraries. TODO: package mgmnt (import, +namespace)
addpath(genpath('./'))

rootData = ['C:/Users/natha/OneDrive - Universiteit Twente/DecompositionAssignment/Data/BioRob (preliminary study)/Subject0_test Michele/'];

%% Load EMG training data and preprocess it 
rawdataPath = 'RT-HDEMG/second/emgs/walking3 - 20231215T153410/';
files = {'walking3-20231215T153410.DATA.Poly5'};

traindata_obj= TMSiSAGA.Poly5.read([rootData,rawdataPath,files{1}]);

fs=    traindata_obj.sample_rate;
ch_names =   {cell2mat(traindata_obj.channels).alternative_name}';
monopolarIDs = contains({cell2mat(traindata_obj.channels).alternative_name},'G');
ch_names = ch_names(monopolarIDs);
monopolar_EMGs  =traindata_obj.samples(monopolarIDs,:);
syncChannel1 =traindata_obj.samples(67,:);  % 60s (13782 to 136650) -200  (13582 to 136450)   23582 123432
syncChannel = syncChannel1(13782:136651);
monopolar_EMGs  = monopolar_EMGs(:,13782:136651);

EMGtrain = monopolar_EMGs;
% discards noisy channels
[~, ~, ~,emgMask] = EMGfilter(EMGtrain, 'discardChannels','var','fcNotch' ,[]); % 'fcNotch' ,[]
% substracts average referencing
EMGtrain = re_referencing(EMGtrain,'referenceConfig' ,{1:32},'emgMask',emgMask,'plotFlag',0);
% filtering (band-pass, notch)
[EMGtrain, EMGenvelope, filterParam,~,qualityResults] = EMGfilter(EMGtrain,'checkNotch',true); %'channelSelection',ch); % 'fcNotch' ,[]

%savePath = 'EMGs\Processed\';
%save([savePath,'walking_filtered.mat'],'EMGtrain','emgMask','filterParam')

clear monopolar_EMGs testdata_obj

%% Offline decomposition training data
% 10 iterations to keep the process fast. No tanh denoising is used as this
% is only applicable on whitened data. The results filters might be a
% little bit biased towards the higher peaks in the signal. 

for i=1%:length(files)
    %load([loadPath,files{i}])
    [signaltrain,decompParameters_train] = decompFastICA_v2(EMGtrain,...
        'fs',fs,'emgMask',emgMask,'nbIterations', 10, ...
        'nbextchan', 1000,'qc_threshold', 0.85,'showPlots',1); %'tanh_denoise', 3
    %'tanh_denoise',3
    %signaltrain.timeStamps =  timeStamps;
    %save([savePath,files{i}],'signal','decompParameters')
end

%% Extract online parameters of the training data
% whiten = 0 equals the CKC algorithm with non whitened data. 
% tanh_denoise = 3 means that if whiten = 1, tanh denoising is done on the
% whitened extended data with std = 3;
% inputParameters function extracts the correct (non) whitened parameters
% depending on the input from whiten and tanh
% cut = 1 because in decompFastICA_v21 the signal is also cut
whiten = 0;
tanh_denoise = 0;
onlineDecomp = true;
cut = 1;
if onlineDecomp
    fs = 2048;
    [eSIGtrain, parameters_train] = inputParameters(EMGtrain, decompParameters_train, signaltrain, whiten, tanh_denoise, fs,cut);
    parameters_train.buffersize = 100;
    [OnlineDecompParameters_train] = getonlineparameters(eSIGtrain, parameters_train, fs);
    [OnlineDecompParameters_train.PulseT,OnlineDecompParameters_train.Distime, ~] = getPulseT(eSIGtrain, parameters_train, fs);
    OnlineDecompParameters_train.EMGmask = decompParameters_train.EMGmask;
end 

%% Get common discharges between whitened and non-whitened calculation of PulseT
[comdis,roa] = getcomdis(signaltrain.PulseT, OnlineDecompParameters_train.PulseT, signaltrain.Distime , OnlineDecompParameters_train.Distime, round(fs/40), 0.00025, 0.3, fs);
%% Compare non whitened and whitened spike trains of training data.
close all
plot(OnlineDecompParameters_train.PulseT(1,:));
hold on
plot(signaltrain.PulseT(1,:));
hold off
legend('Non Whitened', 'whitened');

%% Load EMG test data and preprocess it 
rawdataPath = 'RT-HDEMG/second/emgs/walking1 - 20231215T152949/';
files = {'walking1-20231215T152949.DATA.Poly5'};

testdata_obj= TMSiSAGA.Poly5.read([rootData,rawdataPath,files{1}]);

fs=    testdata_obj.sample_rate;
ch_names =   {cell2mat(testdata_obj.channels).alternative_name}';
monopolarIDs = contains({cell2mat(testdata_obj.channels).alternative_name},'G');
ch_names = ch_names(monopolarIDs);
monopolar_EMGs  =testdata_obj.samples(monopolarIDs,:);
syncChannel =testdata_obj.samples(67,:);  % 60s (13782 to 136650) -200  (13582 to 136450)   23582 123432
syncChannel = syncChannel(13782:136651);
monopolar_EMGs  = monopolar_EMGs(:,13782:136651);

EMGtest = monopolar_EMGs;
% discards noisy channels
%[~, ~, ~,emgMask] = EMGfilter(EMGtest, 'discardChannels','var','fcNotch' ,[]); % 'fcNotch' ,[]
% substracts average referencing

EMGtest = re_referencing(EMGtest,'referenceConfig' ,{1:32},'emgMask',emgMask,'plotFlag',0);
% filtering (band-pass, notch)
[EMGtest, EMGenvelope1, filterParam1,~,qualityResults1] = EMGfilter(EMGtest,'checkNotch',true); %'channelSelection',ch); % 'fcNotch' ,[]

%savePath = 'EMGs\Processed\';
%save([savePath,'walking_filtered.mat'],'EMG','emgMask','filterParam')

clear monopolar_EMGs testdata_obj

%% Offline decomposition test data
% 10 iterations to keep the process fast. No tanh denoising is used as this
% is only applicable on whitened data. The results filters might be a
% little bit biased towards the higher peaks in the signal. 

for i=1%:length(files)
    [signaltest,decompParameters_test] = decompFastICA_v21(EMGtest,...
        'fs',fs,'emgMask',emgMask,'nbIterations', 10, ...
        'nbextchan', 1000,'SILthreshold', 0.85,'showPlots',1); %'tanh_denoise', 3
    %'tanh_denoise',3
    %signaltest.timeStamps =  timeStamps;
    %save([savePath,files{i}],'signal','decompParameters')
end

%% Extract online parameters of the test data
% whiten = 0 equals the CKC algorithm with non whitened data. 
% tanh_true = 1 means that if whiten = 1, tanh denoising is done on the whitened extended data. 
% inputParameters function extracts the correct (non) whitened parameters
% depending on the input from whiten and tanh

whiten = 0;
tanh_denoise = 0;
onlineDecomp = true;
cut = 1;
if onlineDecomp
    fs = 2048;
    [eSIGtest, parameters_test] = inputParameters(EMGtest, decompParameters_test, signaltest, whiten, tanh_denoise, fs,cut);
    parameters_test.buffersize = 100;
    [OnlineDecompParameters_test] = getonlineparameters(eSIGtest, parameters_test, fs);
    [OnlineDecompParameters_test.PulseT,OnlineDecompParameters_test.Distime, ~] = getPulseT(eSIGtest, parameters_test, fs);
    OnlineDecompParameters_test.EMGmask = decompParameters_test.EMGmask;
end 

%% Psuedo online decomposition
% Uses the found filters from EMGtrain data on EMGtest data 
[signalAdapt,AdaptDecompParameters] = pseudoAdaptDecomp(EMGtest, OnlineDecompParameters_train, 'refreshRate', 1); %'tanh_denoise', 3
load([rootData, 'train_test_data/adapt_data.mat']);

%% Get common discharges and Rate of Agreement of offline spike trains and online adaptive spike trains.
[comdis,roa] = getcomdis(OnlineDecompParameters_test.PulseT(:, 1:size(signalAdapt.PulseT,2)), signalAdapt.PulseT, OnlineDecompParameters_test.Distime , signalAdapt.Distime, round(fs/40), 0.00025, 0.3, fs);
%% Compare outcomes of offline decomposition and online adaptive decomposition
close all
plot(signalAdapt.PulseT(1,:));
hold on
plot(OnlineDecompParameters_test.PulseT(1, :));
%plot(EMGtest(20,fs:end-fs)/mean(maxk(EMGtest(20,fs:end-fs),10)));
hold off
legend('Adaptive', 'non-adaptive')
%%
getFiringProperties({signalAdapt.Distime{1},signaltest.Distime{1}},'flagPlotDR',1);
legend('Adapt','Non Adapt');

%% Save decomposition results
save([rootData, 'train_test_data/train_data.mat'], 'OnlineDecompParameters_train', 'signaltrain', 'decompParameters_train');
save([rootData, 'train_test_data/test_data.mat'], 'OnlineDecompParameters_test', 'signaltest', 'decompParameters_test');
save([rootData, 'train_test_data/adapt_data.mat'], 'AdaptDecompParameters', "signalAdapt");


