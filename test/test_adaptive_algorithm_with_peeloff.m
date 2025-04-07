%% Same as test_adaptive_algorithm except for some extra peel off steps
% in both test and training sets, this will produce more MUfilters, as 1
% peel off step is implemented. 
s
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

rootData = ['C:/Users/natha/OneDrive - Universiteit Twente/DecompositionAssignment/Data/BioRob (preliminary study)/Subject0_test Michele/'];
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


clear monopolar_EMGs testdata_obj

%% Offline decomposition training data
% 10 iterations to keep the process fast. No tanh denoising is used as this
% is only applicable on whitened data. The results filters might be a
% little bit biased towards the higher peaks in the signal. 

for i=1%:length(files)
    %load([loadPath,files{i}])
    [signaltrain,decompParameters_train] = decompFastICA_v21(EMGtrain,...
        'fs',fs,'emgMask',emgMask,'nbIterations', 20, ...
        'nbextchan', 1000,'SILthreshold', 0.85,'showPlots',1); %'tanh_denoise', 3
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

whiten = 0;
tanh_denoise = 0;
onlineDecomp = true;
cut = 1;
if onlineDecomp
    fs = 2048;
    [eSIGtrain, parameters_train] = inputParameters(EMGtrain, decompParameters_train, signaltrain, whiten, tanh_denoise, fs, cut);
    parameters_train.buffersize = 50;
    [OnlineDecompParameters_train] = getonlineparameters(eSIGtrain, parameters_train, fs);
    [OnlineDecompParameters_train.PulseT,OnlineDecompParameters_train.Distime, ~] = getPulseT(eSIGtrain, parameters_train, fs);
    OnlineDecompParameters_train.EMGmask = decompParameters_train.EMGmask;
end 
%% Peel off training data (optional)
muapWin = round(.0255*fs); %ms

[EMGtrainres,What_train]=...
    peeloff_LS(signaltrain.spikeTrains(1:length(EMGtrain(:,fs:end-fs)),:), EMGtrain(:,fs:end-fs), muapWin);
%% Offline decomposition training Peel-off data (optional)
% 10 iterations to keep the process fast. No tanh denoising is used as this
% is only applicable on whitened data. The results filters might be a
% little bit biased towards the higher peaks in the signal. 

for i=1%:length(files)
    %load([loadPath,files{i}])
    [signaltrain_res,decompParameters_train_res] = decompFastICA_v21(EMGtrainres,...
        'fs',fs,'emgMask',emgMask,'nbIterations', 5, ...
        'nbextchan', 1000,'SILthreshold', 0.85,'showPlots',1, 'cut', 0); %'tanh_denoise', 3
    %'tanh_denoise',3
    %signaltrain.timeStamps =  timeStamps;
    %save([savePath,files{i}],'signal','decompParameters')
end

% Example parallelDeocmposition (useful for more than 4 files)
% pendingfiles  = parallelDecomposition(files, loadPath, savePath, 'fs',2048,...
%     'filterEMGsflag',0,'nbextchan',1000,'SILthreshold', 0.85,'nbIterations', 50,...
%     'emgMask',{emgMask}); %  'timeStamps',timeStamps
%% Extract online parameters of the peel off training data (Optional)
% whiten = 0 equals the CKC algorithm with non whitened data. 
% tanh_denoise = 3 means that if whiten = 1, tanh denoising is done on the
% whitened extended data with std = 3;
% inputParameters function extracts the correct (non) whitened parameters
% depending on the input from whiten and tanh

whiten = 0;
tanh_denoise = 0;
onlineDecomp = true;
cut = 0;
if onlineDecomp
    fs = 2048;
    [eSIGtrain_res, parameters_train_res] = inputParameters(EMGtrainres, decompParameters_train_res, signaltrain_res, whiten, tanh_denoise, fs, cut);
    parameters_train_res.buffersize = 50;
    [OnlineDecompParameters_train_res] = getonlineparameters(eSIGtrain_res, parameters_train_res, fs);
    [OnlineDecompParameters_train_res.PulseT,OnlineDecompParameters_train_res.Distime, ~] = getPulseT(eSIGtrain_res, parameters_train_res, fs);
    OnlineDecompParameters_train_res.EMGmask = decompParameters_train.EMGmask;
end 


%% Get common discharges between whitened and non-whitened calculation of PulseT
[comdis] = getcomdis(signaltrain_comb.PulseT, OnlineDecompParameters_train_comb.PulseT, signaltrain_comb.Distime , OnlineDecompParameters_train_comb.Distime, round(fs/40), 0.00025, 0.3, fs);
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
        'fs',fs,'emgMask',emgMask,'nbIterations', 50, ...
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
if onlineDecomp
    fs = 2048;
    [eSIGtest, parameters_test] = inputParameters(EMGtest, decompParameters_test, signaltest, whiten, tanh_denoise, fs);
    parameters_test.buffersize = 50;
    [OnlineDecompParameters_test] = getonlineparameters(eSIGtest, parameters_test, fs);
    [OnlineDecompParameters_test.PulseT,OnlineDecompParameters_test.Distime, ~] = getPulseT(eSIGtest, parameters_test, fs);
    OnlineDecompParameters_test.EMGmask = decompParameters_test.EMGmask;
end 

%% Peel off test data (optional)
muapWin = round(.0255*fs); %ms

[EMGtestres,What_test]=...
    peeloff_LS(signaltest.spikeTrains(1:length(EMGtest(:,fs:end-fs)),:), EMGtest(:,fs:end-fs), muapWin);
%% Offline decomposition training Peel-off data (optional)
% 10 iterations to keep the process fast. No tanh denoising is used as this
% is only applicable on whitened data. The results filters might be a
% little bit biased towards the higher peaks in the signal. 

for i=1%:length(files)
    %load([loadPath,files{i}])
    [signaltest_res,decompParameters_test_res] = decompFastICA_v21(EMGtestres,...
        'fs',fs,'emgMask',emgMask,'nbIterations', 10, ...
        'nbextchan', 1000,'SILthreshold', 0.85,'showPlots',1, 'cut', 0); %'tanh_denoise', 3
    %'tanh_denoise',3
    %signaltrain.timeStamps =  timeStamps;
    %save([savePath,files{i}],'signal','decompParameters')
end

%% Extract online parameters of the peel off training data (Optional)
% whiten = 0 equals the CKC algorithm with non whitened data. 
% tanh_denoise = 3 means that if whiten = 1, tanh denoising is done on the
% whitened extended data with std = 3;
% inputParameters function extracts the correct (non) whitened parameters
% depending on the input from whiten and tanh

whiten = 0;
tanh_denoise = 0;
onlineDecomp = true;
cut = 0;
if onlineDecomp
    fs = 2048;
    [eSIGtest_res, parameters_test_res] = inputParameters(EMGtestres, decompParameters_test_res, signaltest_res, whiten, tanh_denoise, fs, cut);
    parameters_test_res.buffersize = 50;
    [OnlineDecompParameters_test_res] = getonlineparameters(eSIGtest_res, parameters_test_res, fs);
    [OnlineDecompParameters_test_res.PulseT,OnlineDecompParameters_test_res.Distime, ~] = getPulseT(eSIGtest_res, parameters_test_res, fs);
    OnlineDecompParameters_test_res.EMGmask = decompParameters_train.EMGmask;
end 

%% Psuedo online decomposition
% Uses the found filters from EMGtrain data on EMGtest data 
% First apply normal MU filters
[signalAdapt,AdaptDecompParameters] = pseudoAdaptDecomp(EMGtest, OnlineDecompParameters_train, 'refreshRate', 1); %'tanh_denoise', 3
% Apply peel off data trained MU filters
[signalAdapt_res,AdaptDecompParameters_] = pseudoAdaptDecomp(EMGtest, OnlineDecompParameters_train_res, 'refreshRate', 1); %'tanh_denoise', 3

%% Merge structs for combined results
signalAdapt_comb = merge_structs(signalAdapt, signalAdapt_res);
OnlineDecompParameters_test_comb = merge_structs(OnlineDecompParameters_test, OnlineDecompParameters_test_res);
%% Get common discharges and Rate of Agreement of offline spike trains and online adaptive spike trains.
% In this case, the columns represent the ith pulse train of adaptive algorithm
% the rows present the jth pulse train of the offline decomposition results.
[comdis,roa] = getcomdis(OnlineDecompParameters_test_comb.PulseT(:, 1:size(signalAdapt_comb.PulseT,2)), signalAdapt_comb.PulseT, OnlineDecompParameters_test_comb.Distime , signalAdapt_comb.Distime, round(fs/40), 0.00025, 0.3, fs);
%% Compare outcomes of offline decomposition and online adaptive decomposition
close all
% Plots the 3 th pulse train of adaptive algorithm
plot(signalAdapt_comb.PulseT(6,:));
hold on
% Plots the 1st pulse train of the offline decomposition
plot(OnlineDecompParameters_test_comb.PulseT(7, :));
%plot(EMGtest(20,fs:end-fs)/mean(maxk(EMGtest(20,fs:end-fs),10)));
hold off
legend('Adaptive', 'non-adaptive')
%%
getFiringProperties({signalAdapt_comb.Distime{6},OnlineDecompParameters_test_comb.Distime{7}},'flagPlotDR',1);
legend('Adapt','Non Adapt');


