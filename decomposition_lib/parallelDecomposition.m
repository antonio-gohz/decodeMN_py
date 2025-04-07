function pendingfiles  =...
            parallelDecomposition(fileNames, filePath, savePath, varargin)



% Default values for optional parameters for EMG filtering 
fs = 2048;
filterEMGsflag = 0;
emgFiltparams = [];  % files with filterParams and emgMasks
flagTimeStamps = 0;
timeStamps=zeros(length(fileNames),2);
        defaultEMGfiltFlag = 0;

% decomp params
nbIterations = 50;
nbextchan = 1000;
SILthreshold = 0.85;
preOptFilters = [];

% 
peelOffFlag = 0;
muapWin = .0255; %ms
emgMask= cell(size(fileNames));

% Process optional input parameters using a loop and switch case
for j = 1:2:length(varargin)
    param = varargin{j};
    value = varargin{j+1};
    switch param
        case 'fs'
            fs = value;
        case 'filterEMGsflag'           %%% emg params
            filterEMGsflag = value;
        case 'emgFiltparams'          
            emgFiltparams = value;
        case 'emgMask'
            emgMask = value;
        case 'timeStamps'          
            timeStamps = value;
        case 'nbIterations'    %%%%%% decomp params
            nbIterations = value;
        case 'nbextchan'
            nbextchan = value;
        case 'SILthreshold'
            SILthreshold = value;
        case 'preOptFilters'
            preOptFilters = value;
        case 'peelOffFlag'
            peelOffFlag = value;
        case 'muapWin'
            muapWin = value;
        otherwise
            error('Invalid optional input parameter: %s', param);
    end
end

N= length(fileNames);
p=1;

% todo develop function get filter params
if filterEMGsflag  % if emg filter
    if isempty(emgFiltparams)  % if filer with specific params 
        defaultEMGfiltFlag = 1;
        emgFiltparams= cell(size(fileNames));
    elseif length(emgFiltparams)==1
        emgFiltparams = repelem(emgFiltparams,N);
    end
end

if sum(timeStamps)>0 % if  timesptamps
    flagTimeStamps=1;
    if size(timeStamps,2)>size(timeStamps,1) && numel(timeStamps)>2
        timeStamps = timeStamps';
    end
    if length(timeStamps)==2 
        timeStamps= repmat(timeStamps,[N,1]);
    end
    %timeStamps = round(timeStamps*fs);
end


%save([savePath,'timeStampsDecomp'],'timeStamps','fileNames')

muapWin = round(muapWin*fs);

pp = parallel.pool.DataQueue;
h = waitbar(0, 'Please wait ...');
afterEach(pp, @nUpdateWaitbar);

errorVec = zeros(size(fileNames));
numMUs = zeros(size(fileNames));
mean_sil =zeros(size(fileNames));
se_sil = zeros(size(fileNames));
peeloffResults =cell(size(fileNames));
signal=cell(size(fileNames));
decompParameters=cell(size(fileNames));
EMGs=cell(size(fileNames));
struct_name=cell(size(fileNames));
%CE = cell(size(fileNames));

nWorkers = min([N+1,12]);
parfor (i=1:N,nWorkers)
%for i=1:length(fileNames) 
    try
        %load EMGs
        EMGs{i} = load([char(filePath),fileNames{i}],'EMG');

        if flagTimeStamps
            EMGs{i}.EMG = EMGs{i}.EMG(:,timeStamps(i,1):timeStamps(i,2));
        end

        if filterEMGsflag
            if defaultEMGfiltFlag
                EMGs{i}.EMG = EMGfilter(EMGs{i}.EMG, 'fs',fs);
            else
                EMGs{i}.EMG = EMGfilter(EMGs{i}.EMG, 'fs',emgFiltparams{i}.filterParams.fs,'fcNotch',emgFiltparams{i}.filterParams.fcNotch,...
                    'fcBP',emgFiltparams{i}.filterParams.fcBP,'fcEnv',emgFiltparams{i}.filterParams.fcEnv,...
                    'orderFilt',emgFiltparams{i}.filterParams.orderFilt,'q',emgFiltparams{i}.filterParams.q);
            end
        end
        if isempty(emgMask{i})
            emgMask{i}= true(1,min(size(EMGs{i}.EMG)));
        end
        [signal{i},decompParameters{i}] = decompFastICA(EMGs{i}.EMG,'fs',fs,...
            'emgMask',emgMask{i},'nbIterations', nbIterations, ...
            'nbextchan', nbextchan,'SILthreshold', SILthreshold,...
            'preOptFilters',preOptFilters);

        signal{i}.spikeTrains = ...
            getFiringProperties(signal{i}.Dischargetimes,...
            'tVec',linspace(0,length(EMGs{i}.EMG)/fs,length(EMGs{i}.EMG)));
        % save
        struct_name{i} = struct('signal',signal(i),'decompParameters',decompParameters(i));
        parsave([savePath,erase(fileNames{i},'.mat'),'_decomposed.mat'],struct_name{i})
        
        if peelOffFlag
            [peeloffResults{i}.EMGres,peeloffResults{i}.What]=...
                peeloff_LS(signal{i}.spikeTrains, EMGs{i}.EMG, muapWin);
            % save
            struct_name{i} = struct('peeloffResults',peeloffResults{i});
            parsave([savePath,erase(fileNames{i},'.mat'),'_EMGres.mat'],struct_name{i})
        end
        
        numMUs(i) = length(signal{i}.SIL);
        mean_sil(i) = nanmean(signal{i}.SIL);
        se_sil(i) = nanstd(signal{i}.SIL)/sqrt(numMUs(i));

        % clear variables
        signal{i} = [];
        decompParameters{i} = [];
        struct_name{i} =[];
        peeloffResults{i} = [];
        EMGs{i} = [];
        emgFiltparams{i} =[];
    catch 
        send(pp, i+100);
        errorVec(i) = 1;
        % clear variables
        signal{i} = [];
        decompParameters{i} = [];
        struct_name{i} =[];
        peeloffResults{i} = [];
        EMGs{i} = [];
        emgFiltparams{i} =[];
    end
            send(pp, i);    

end
try
pendingfiles = fileNames(logical(errorVec));
plotSummaryQualityMetrics([numMUs',mean_sil'], {'# MUs','Mean SIL'},erase(fileNames,'.mat'), 'metricThresholds', [5,0.85])
savefig(gcf,[savePath,'\qualityResults_decomposition'])
catch
end
    function nUpdateWaitbar(~)
        waitbar(p/N, h);
        p = p + 1;
    end
end