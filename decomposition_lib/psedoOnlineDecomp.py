# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m

    
@function
def psedoOnlineDecomp(EMGs=None,decompParameters=None,varargin=None,*args,**kwargs):
    varargin = psedoOnlineDecomp.varargin
    nargin = psedoOnlineDecomp.nargin

    # simulatedRTdecomp performs real-time decomposition of EMG signals.
# It returns decomposed signals, spike trains, and activation dynamics.
    
    # Inputs:
#   EMGs: EMG signals (columns represent muscles, rows represent samples)
#   decompParameters: Decomposition parameters
#   varargin: Optional input parameters
#       - 'fs': Sampling frequency (default: 2048)
#       - 'idxmuscle1': Index of the muscle to analyze (default: 1)
#       - 'refreshRate': Refresh rate (default: 64) [Hz]
#       - 'flagVis': Visualization flag (default: 0)
#       - 'SmoothingEditFieldDR': Smoothing factor for DR (default: 0.2)[s]
#       - 'plotTimeRange': Time range for plotting (default: 10)[s]
#       - 'flagAD': Activation dynamics flag (default: 0)
#       - 'aps': Peak amplitudes for Activation dynamics 
#       - 'tcs': contaction times Activation dynamics [s]
#       - 'stateMask4MUFilters':
    
    # Outputs:
#   signal: Decomposed EMG signals and spike trains
#   decompParameters: Decomposition parameters
#   activation: Activation dynamics
    
    # Check if EMGs needs to be transposed
    if size(EMGs,1) > size(EMGs,2):
        EMGs=EMGs.T
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:27
    
    # default parameters
    fs=2048
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:31
    idxmuscle1=1
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:32
    refreshRate=10
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:33
    
    decomp_method=2
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:34
    flagAD=0
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:35
    aps=[]
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:36
    tcs=[]
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:37
    flagVis=0
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:39
    SmoothingEditFieldDR=0.5
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:40
    
    plotTimeRange=10
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:41
    
    axesLineConfig=concat([[1,1],[1,2],[2,3]])
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:42
    stateMask4MUFilters=[]
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:45
    for i in arange(1,length(varargin),2).reshape(-1):
        param=varargin[i]
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:48
        value=varargin[i + 1]
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:49
        if 'fs' == param:
            fs=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:52
        else:
            if 'idxmuscle1' == param:
                idxmuscle1=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:54
            else:
                if 'refreshRate' == param:
                    refreshRate=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:56
                else:
                    if 'flagAD' == param:
                        flagAD=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:58
                    else:
                        if 'flagVis' == param:
                            flagVis=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:60
                        else:
                            if 'decomp_method' == param:
                                decomp_method=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:62
                            else:
                                if 'SmoothingEditFieldDR' == param:
                                    SmoothingEditFieldDR=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:64
                                else:
                                    if 'plotTimeRange' == param:
                                        plotTimeRange=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:66
                                    else:
                                        if 'aps' == param:
                                            aps=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:68
                                        else:
                                            if 'tcs' == param:
                                                tcs=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:70
                                            else:
                                                if 'stateMask4MUFilters' == param:
                                                    stateMask4MUFilters=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:72
                                                else:
                                                    if 'tanhNormEMG' == param:
                                                        tanhNormEMG=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:74
                                                    else:
                                                        if 'tanhNormWEMG' == param:
                                                            tanhNormWEMG=copy(value)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:76
                                                        else:
                                                            error('Invalid optional input parameter: %s',param)
    
    if logical_not(isempty(stateMask4MUFilters)):
        decompParameters.MUfiltxiReSIGt[arange(),logical_not(stateMask4MUFilters)]=[]
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:83
        decompParameters.MUFilters[arange(),logical_not(stateMask4MUFilters)]=[]
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:84
        decompParameters.norm[logical_not(stateMask4MUFilters)]=[]
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:85
        decompParameters.centroid[logical_not(stateMask4MUFilters),arange()]=[]
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:87
        MUfilt=decompParameters.MUfiltxiReSIGt.T
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:88
    else:
        MUfilt=dot(decompParameters.MUFilters.T,decompParameters.iReSIGt)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:90
    
    # Simulating adquisition parameters
    ComParameters.nsamp = copy(floor(fs / refreshRate))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:93
    
    OnlineParameters.durationonline = copy(length(EMGs))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:94
    nwin=floor(OnlineParameters.durationonline / ComParameters.nsamp)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:95
    
    # Pre allocate empty matrices for decomposition
    EMG=zeros(size(EMGs))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:98
    time_DR=linspace(0,size(EMG,2) / fs,nwin)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:99
    time_spikes_AD=linspace(0,size(EMG,2) / fs,size(EMG,2))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:100
    EMGtmp=zeros(size(EMG,1),ComParameters.nsamp)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:101
    esample2=zeros(dot((size(EMG,1) - sum(logical_not(decompParameters.EMGmask))),decompParameters.extensionfactor(idxmuscle1)),decompParameters.extensionfactor(idxmuscle1) - 1)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:102
    signal.spikeTrains = copy(zeros(size(decompParameters.MUFilters,2),dot(nwin,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:103
    signal.Pulsetrain = copy(zeros(size(decompParameters.MUFilters,2),dot(ComParameters.nsamp,nwin)))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:104
    noise_centroids=multiply(ones(ComParameters.nsamp,size(decompParameters.MUFilters,2)),(decompParameters.centroid(arange(),1).T))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:105
    spike_centroids=multiply(ones(ComParameters.nsamp,size(decompParameters.MUFilters,2)),(decompParameters.centroid(arange(),2).T))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:106
    #MUfilt = decompParameters.MUFilters' * decompParameters.iReSIGt;
    
    DR=zeros(nwin,size(decompParameters.MUFilters,2))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:111
    ct_dec=zeros(concat([1,nwin]))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:112
    # Pre allocate empty matrices for activation dynamics (AD)
    muActivation=NaN(size(EMG,2),size(decompParameters.MUFilters,2))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:115
    totalActivation=NaN(size(EMG,1),1)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:116
    ct_AD=zeros(concat([1,nwin]))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:117
    intialValues=zeros(concat([2,size(decompParameters.MUFilters,2)]))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:118
    if flagAD:
        axesLineConfig=concat([[1,1],[1,2],[2,3],[3,4]])
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:121
    
    if flagVis:
        SmoothingEditFieldDR=round(dot(SmoothingEditFieldDR,refreshRate))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:128
        plotLengthWindows=round(dot(refreshRate,plotTimeRange) / 2)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:129
        numMUs=size(decompParameters.MUFilters,2)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:130
        h,ax,__=preparePlots(numMUs,'axesLineConfig',axesLineConfig,nargout=3)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:131
    
    k=1
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:134
    f=waitbar(0,'Decomposition...')
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:135
    while k <= nwin:

        EMG[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=EMGs(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:137
        tic
        if decomp_method == 1:
            EMGtmp=EMG(decompParameters.EMGmask,arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:141
            signal.Pulsetrain(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))),signal.spikeTrains(arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)),arange()),esample2=getspikesonline(EMGtmp,decompParameters.extensionfactor(idxmuscle1),esample2,decompParameters.MUFilters,decompParameters.whitMat,decompParameters.normIPT,decompParameters.centroid,noise_centroids,spike_centroids,ComParameters.nsamp,fs,decompParameters.normtanh,nargout=3)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:142
        else:
            if decomp_method == 2:
                EMGtmp=EMG(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:146
                signalTmp,decompParameters=getspikes_online(EMGtmp,decompParameters,fs,true,nargout=2)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:147
                signal.Pulsetrain[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.Pulsetrain
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:148
                signal.spikeTrains[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.spikeTrains
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:149
            else:
                if decomp_method == 3:
                    EMGtmp=EMG(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:151
                    X=extend(EMGtmp(decompParameters.EMGmask,arange()),decompParameters.extensionfactor)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:152
                    X=X(arange(),arange(1,size(EMGtmp,2)))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:153
                    decompParameters.norm = copy(decompParameters.normIPT)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:154
                    signal.Pulsetrain(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))),signal.spikeTrains(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))),decompParameters,esample2=getspikesonline_Yeung(X,decompParameters,[],fs,nargout=4)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:157
        ct_dec[k]=toc
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:160
        if flagAD:
            tic
            muActivation(arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)),arange()),totalActivation(arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))=muActivationFuglevand(signal.spikeTrains(arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)),arange()),aps,tcs,intialValues,nargout=2)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:165
            intialValues=muActivation(arange(dot((k),ComParameters.nsamp) - 1,dot((k),ComParameters.nsamp)),arange())
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:168
            ct_AD[k]=toc
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:169
        if flagVis:
            if k > SmoothingEditFieldDR:
                DR[k,arange()]=dot(sum(signal.spikeTrains(arange(dot((k - SmoothingEditFieldDR),ComParameters.nsamp),dot(k,ComParameters.nsamp)),arange())),fs) / (dot(ComParameters.nsamp,SmoothingEditFieldDR))
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:174
                DR[k,arange()]=mean(DR(arange(k - round(SmoothingEditFieldDR / 3),k),arange()),1)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:175
            if k < plotLengthWindows:
                updatePlot(cellarray([[time_DR(k),DR(k,arange())],[time_DR(arange(1,k)),DR(arange(1,k),arange())],[time_spikes_AD(arange(1,dot(k,ComParameters.nsamp))),(dot(0.8,(signal.spikeTrains(arange(1,dot(k,ComParameters.nsamp)),arange()))) + (arange(1,numMUs)) - 0.4)],[time_spikes_AD(arange(1,dot(k,ComParameters.nsamp))),muActivation(arange(1,dot(k,ComParameters.nsamp)),arange()) / 7.5]]),h,ax,concat([time_DR(1),time_DR(dot(plotLengthWindows,2))]),concat([0,1,1]))
            else:
                updatePlot(cellarray([[time_DR(k),DR(k,arange())],[time_DR(arange(k - plotLengthWindows + 1,k)),DR(arange(k - plotLengthWindows + 1,k),arange())],[time_spikes_AD(arange(dot((k - plotLengthWindows + 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))),(dot(0.8,(signal.spikeTrains(arange(dot((k - plotLengthWindows + 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)),arange()))) + (arange(1,numMUs)) - 0.4)],[time_spikes_AD(arange(dot((k - plotLengthWindows + 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))),muActivation(arange(dot((k - plotLengthWindows + 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)),arange()) / 7.5]]),h,ax,concat([time_DR(k - plotLengthWindows + 1),time_DR(k) + plotTimeRange / 2]),concat([0,1,1]))
            pause(ComParameters.nsamp / fs)
        tit=concat(['Block : ',num2str(k),' of ',num2str(nwin)])
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:194
        waitbar((k) / (nwin),f,tit)
        k=k + 1
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:196

    
    for m in arange(1,size(signal.spikeTrains,2)).reshape(-1):
        signal.Dischargetimes[m]=find(signal.spikeTrains(arange(),m).T > 0)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:200
    
    # Store signal and activation variables
    signal.spikeTrains = copy(signal.spikeTrains)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:204
    signal.ct = copy(ct_dec)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:205
    activation.muActivation = copy(muActivation)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:207
    activation.totalActivation = copy(totalActivation)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:208
    activation.ct = copy(ct_AD)
# ..\decompositionproject_python\decomposition_lib\psedoOnlineDecomp.m:209
    return signal,decompParameters,activation
    
if __name__ == '__main__':
    pass
    