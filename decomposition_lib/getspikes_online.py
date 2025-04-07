# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m

    
@function
def getspikes_online(EMGs=None,decompParameters=None,fs=None,adaptiveFlag=None,comp_method=None,esample2=None,*args,**kwargs):
    varargin = getspikes_online.varargin
    nargin = getspikes_online.nargin

    # simulatedRTdecomp performs real-time decomposition of EMG signals.
# It returns decomposed signals, spike trains, and activation dynamics.
    
    # Inputs:
#   EMGs: EMG signals (columns represent muscles, rows represent samples)
#   decompParameters: Decomposition parameters
#   varargin: Optional input parameters
#       - 'fs': Sampling frequency (default: 2048)
    
    # Outputs:
#   signal: Decomposed EMG signals and spike trains
#   decompParameters: Decomposition parameters
    
    # Check if EMGs needs to be transposed
# if size(EMGs, 1) > size(EMGs, 2)
#     EMGs = EMGs';
# end
    
    #
    if nargin < 4:
        adaptiveFlag=copy(false)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:22
    
    if nargin < 5:
        comp_method=1
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:25
    
    tanh_denoise=3
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:28
    lambda_=0.985
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:29
    zvalue=3.3
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:30
    nsamp=size(EMGs,2)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:31
    esample1=extend(EMGs(decompParameters.EMGmask,arange()),decompParameters.extensionfactor)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:33
    X=esample1(arange(),arange(1,nsamp))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:34
    X[arange(),arange(1,decompParameters.extensionfactor - 1)]=X(arange(),arange(1,decompParameters.extensionfactor - 1)) + esample2
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:35
    esample2=esample1(arange(),arange(nsamp + 1,end()))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:36
    # 
# esample1 = extend(EMGtmp, extensionfactor);
#     esample = esample1(:,1:nsamp);
#     esample(:,1:extensionfactor-1) = esample(:,1:extensionfactor-1) + esample2;
#     esample2 = esample1(:,nsamp+1:end);
    
    # whitened MU filter * whitened data
    
    #X = decompParameters.whitMat*eSIG;
#normtanh =  tanh_denoise*std(X,[],2);
#if ~isempty(decompParameters.normtanh)
    #X = decompParameters.normtanh.*tanh(X./decompParameters.normtanh);
#    X = normtanh.*tanh(X./normtanh);
    
    #end
    if adaptiveFlag:
        ReSIG_win=dot(X,X.T) / length(X)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:54
        decompParameters.ReSIG = copy(dot(decompParameters.ReSIG,lambda_) + dot((1 - lambda_),ReSIG_win))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:55
        decompParameters.iReSIGt = copy(pinv(decompParameters.ReSIG))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:56
        buffer_size=size(decompParameters.buffer_eSIG,2)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:57
    
    PulseT=(dot(dot(decompParameters.MUFilters.T,decompParameters.iReSIGt),X))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:59
    PulseT=multiply(PulseT,abs(PulseT))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:60
    PulseT=PulseT / decompParameters.normIPT.T
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:61
    PulseT=tanh(PulseT)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:62
    mean_noise=zeros(concat([size(decompParameters.MUFilters,2),1]))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:63
    if comp_method == 1:
        spikesTmp,__=islocalmax(PulseT,2,'MinSeparation',round(dot(fs,0.02)),nargout=2)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:66
    else:
        if comp_method == 2 and logical_not(adaptiveFlag):
            spikesTmp=zeros(size(PulseT))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:68
            for i in arange(1,size(decompParameters.MUFilters,2)).reshape(-1):
                __,spikes_ids=findpeaks(PulseT(i,arange()),'MinPeakDistance',round(dot(fs,0.02)),nargout=2)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:70
                spikesTmp[i,spikes_ids]=true
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:71
        else:
            if comp_method == 2 and adaptiveFlag:
                spikes1=false(size(PulseT))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:74
                peaksnoise=false(size(PulseT))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:75
    
    if adaptiveFlag:
        # See if z-score is above treshold.
        PulseZ=zscore(PulseT,0,2)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:81
        if comp_method == 1:
            spikes1=copy(spikesTmp)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:83
            spikes1[PulseZ < zvalue]=0
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:84
            spikes1[PulseT(spikes1) < 0.1]=0
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:85
            peaksnoise=xor(spikes1,spikesTmp)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:86
        for i in arange(1,size(decompParameters.MUFilters,2)).reshape(-1):
            if comp_method == 2:
                __,spikes_ids=findpeaks(PulseT(i,arange()),'MinPeakDistance',round(dot(fs,0.02)),nargout=2)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:90
                #spikesTmp(isoutlier(PulseT(spikes),'percentiles',[0,99]))=[]; # NEW add remove the outliers
            # See if z-score is above treshold.
                PotSpikes=find(PulseZ(i,spikes_ids) > zvalue)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:93
                peaksnoise_ids=setdiff(spikes_ids,spikes_ids(PotSpikes))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:94
                PotSpikes=PotSpikes(PulseT(spikes_ids(PotSpikes)) > 0.1)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:95
                spikes1[i,PotSpikes]=true
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:96
                peaksnoise[i,peaksnoise_ids]=true
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:97
            buffer_spikes=cat(2,decompParameters.buffer_spikes(i,arange()),PulseT(i,spikes1(i,arange())))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:100
            decompParameters.buffer_spikes[i,arange()]=buffer_spikes(arange(end() - buffer_size + 1,end()))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:101
            buffer_eSIG=cat(2,decompParameters.buffer_eSIG(arange(),arange(),i),X(arange(),spikes1(i,arange())))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:102
            decompParameters.buffer_eSIG[arange(),arange(),i]=buffer_eSIG(arange(),arange(end() - buffer_size + 1,end()))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:103
            decompParameters.MUFilters[arange(),i]=mean(decompParameters.buffer_eSIG(arange(),arange(),i),2,'omitnan')
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:104
            mean_noise[i,arange()]=mean(PulseT(i,peaksnoise(i,arange())))
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:105
        #centroid
        prev_centroid_noise=min(decompParameters.centroid,[],2)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:108
        decompParameters.centroid[arange(),1]=mean(decompParameters.buffer_spikes,2,'omitnan')
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:109
        decompParameters.centroid[arange(),2]=dot(lambda_,prev_centroid_noise) + dot((1 - lambda_),mean_noise)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:110
    else:
        spikes1=multiply(PulseT,spikesTmp) >= abs(diff(decompParameters.centroid,1,2)) / 2 + min(decompParameters.centroid,[],2)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:112
    
    signal.spikeTrains = copy(spikes1)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:114
    signal.Pulsetrain = copy(PulseT)
# ..\decompositionproject_python\decomposition_lib\getspikes_online.m:115