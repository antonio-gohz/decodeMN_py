# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m

    
@function
def simulate_online_decomposition(EMGs=None,decompParameters=None,fs=None,varargin=None,*args,**kwargs):
    varargin = simulate_online_decomposition.varargin
    nargin = simulate_online_decomposition.nargin

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
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:27
    
    # default parameters
    idxmuscle1=1
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:31
    refreshRate=10
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:32
    
    decomp_method=2
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:33
    flagAD=0
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:34
    aps=[]
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:35
    tcs=[]
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:36
    plotResults=copy(false)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:37
    flagVis=0
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:38
    SmoothingEditFieldDR=0.5
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:39
    
    plotTimeRange=10
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:40
    
    axesLineConfig=concat([[1,1],[1,2],[2,3]])
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:41
    for i in arange(1,length(varargin),2).reshape(-1):
        param=varargin[i]
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:46
        value=varargin[i + 1]
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:47
        if 'idxmuscle1' == param:
            idxmuscle1=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:50
        else:
            if 'refreshRate' == param:
                refreshRate=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:52
            else:
                if 'flagAD' == param:
                    flagAD=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:54
                else:
                    if 'flagVis' == param:
                        flagVis=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:56
                    else:
                        if 'decomp_method' == param:
                            decomp_method=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:58
                        else:
                            if 'SmoothingEditFieldDR' == param:
                                SmoothingEditFieldDR=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:60
                            else:
                                if 'plotTimeRange' == param:
                                    plotTimeRange=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:62
                                else:
                                    if 'aps' == param:
                                        aps=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:64
                                    else:
                                        if 'tcs' == param:
                                            tcs=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:66
                                        else:
                                            if 'tanhNormEMG' == param:
                                                tanhNormEMG=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:68
                                            else:
                                                if 'plotResults' == param:
                                                    plotResults=copy(value)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:70
                                                else:
                                                    error('Invalid optional input parameter: %s',param)
    
    # Simulating adquisition parameters
    ComParameters.nsamp = copy(floor(fs / refreshRate))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:77
    
    OnlineParameters.durationonline = copy(length(EMGs))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:78
    nwin=floor(OnlineParameters.durationonline / ComParameters.nsamp)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:79
    
    # Pre allocate empty matrices for decomposition
    EMG=zeros(size(EMGs))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:82
    time_DR=linspace(0,size(EMG,2) / fs,nwin)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:83
    time_spikes_AD=linspace(0,size(EMG,2) / fs,size(EMG,2))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:84
    EMGtmp=zeros(size(EMG,1),ComParameters.nsamp)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:85
    numMUs=size(decompParameters.MUFilters,2)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:86
    esample2=zeros(dot(sum(decompParameters.EMGmask),decompParameters.extensionfactor),decompParameters.extensionfactor - 1)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:87
    signal.spikeTrains = copy(zeros(numMUs,dot(nwin,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:88
    signal.Pulsetrain = copy(zeros(numMUs,dot(ComParameters.nsamp,nwin)))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:89
    noise_centroids=multiply(ones(ComParameters.nsamp,numMUs),(decompParameters.centroid(arange(),1).T))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:90
    spike_centroids=multiply(ones(ComParameters.nsamp,numMUs),(decompParameters.centroid(arange(),2).T))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:91
    DR=zeros(nwin,numMUs)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:93
    ct_dec=zeros(concat([1,nwin]))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:94
    
    # Pre allocate empty matrices for activation dynamics (AD)
    muActivation=NaN(size(EMG,2),numMUs)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:97
    totalActivation=NaN(size(EMG,1),1)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:98
    ct_AD=zeros(concat([1,nwin]))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:99
    intialValues=zeros(concat([2,numMUs]))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:100
    if flagAD:
        axesLineConfig=concat([[1,1],[1,2],[2,3],[3,4]])
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:103
    
    if flagVis:
        SmoothingEditFieldDR=round(dot(SmoothingEditFieldDR,refreshRate))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:110
        plotLengthWindows=round(dot(refreshRate,plotTimeRange) / 2)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:111
        h,ax,__=preparePlots(numMUs,'axesLineConfig',axesLineConfig,nargout=3)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:112
    
    buffer_count=zeros(concat([1,numMUs]))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:114
    k=1
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:115
    f=waitbar(0,'Decomposition...')
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:116
    while k <= nwin:

        EMG[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=EMGs(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:118
        tic
        if decomp_method == 1:
            EMGtmp=EMG(decompParameters.EMGmask,arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:122
            signal.Pulsetrain(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))),signal.spikeTrains(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))),esample2=getspikesonline(EMGtmp,decompParameters.extensionfactor(idxmuscle1),esample2,decompParameters.MUFilters,decompParameters.whitMat,decompParameters.normIPT,decompParameters.centroid,noise_centroids,spike_centroids,ComParameters.nsamp,fs,decompParameters.normtanh,nargout=3)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:123
        else:
            if decomp_method == 2:
                EMGtmp=EMG(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:127
                signalTmp,decompParameters,esample2=getspikes_online(EMGtmp,decompParameters,fs,false,1,esample2,nargout=3)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:128
                signal.Pulsetrain[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.Pulsetrain
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:129
                signal.spikeTrains[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.spikeTrains
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:130
            else:
                if decomp_method == 3:
                    EMGtmp=EMG(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:132
                    signalTmp,decompParameters,esample2=getspikes_online(EMGtmp,decompParameters,fs,false,2,esample2,nargout=3)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:133
                    signal.Pulsetrain[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.Pulsetrain
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:134
                    signal.spikeTrains[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.spikeTrains
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:135
                else:
                    if decomp_method == 4:
                        EMGtmp=EMG(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:137
                        signalTmp,decompParameters,buffer_count=getspikes_online(EMGtmp,decompParameters,fs,true,1,buffer_count,nargout=3)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:138
                        signal.Pulsetrain[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.Pulsetrain
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:139
                        signal.spikeTrains[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.spikeTrains
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:140
                    else:
                        if decomp_method == 5:
                            EMGtmp=EMG(arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:142
                            signalTmp,decompParameters,esample2=getspikes_online(EMGtmp,decompParameters,fs,true,2,esample2,nargout=3)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:143
                            signal.Pulsetrain[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.Pulsetrain
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:144
                            signal.spikeTrains[arange(),arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))]=signalTmp.spikeTrains
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:145
        ct_dec[k]=toc
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:148
        if flagAD:
            tic
            muActivation(arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)),arange()),totalActivation(arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)))=muActivationFuglevand(signal.spikeTrains(arange(dot((k - 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)),arange()),aps,tcs,intialValues,nargout=2)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:153
            intialValues=muActivation(arange(dot((k),ComParameters.nsamp) - 1,dot((k),ComParameters.nsamp)),arange())
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:156
            ct_AD[k]=toc
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:157
        if flagVis:
            if k > SmoothingEditFieldDR:
                DR[k,arange()]=dot(sum(signal.spikeTrains(arange(dot((k - SmoothingEditFieldDR),ComParameters.nsamp),dot(k,ComParameters.nsamp)),arange())),fs) / (dot(ComParameters.nsamp,SmoothingEditFieldDR))
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:162
                DR[k,arange()]=mean(DR(arange(k - round(SmoothingEditFieldDR / 3),k),arange()),1)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:163
            if k < plotLengthWindows:
                updatePlot(cellarray([[time_DR(k),DR(k,arange())],[time_DR(arange(1,k)),DR(arange(1,k),arange())],[time_spikes_AD(arange(1,dot(k,ComParameters.nsamp))),(dot(0.8,(signal.spikeTrains(arange(1,dot(k,ComParameters.nsamp)),arange()))) + (arange(1,numMUs)) - 0.4)],[time_spikes_AD(arange(1,dot(k,ComParameters.nsamp))),muActivation(arange(1,dot(k,ComParameters.nsamp)),arange()) / 7.5]]),h,ax,concat([time_DR(1),time_DR(dot(plotLengthWindows,2))]),concat([0,1,1]))
            else:
                updatePlot(cellarray([[time_DR(k),DR(k,arange())],[time_DR(arange(k - plotLengthWindows + 1,k)),DR(arange(k - plotLengthWindows + 1,k),arange())],[time_spikes_AD(arange(dot((k - plotLengthWindows + 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))),(dot(0.8,(signal.spikeTrains(arange(dot((k - plotLengthWindows + 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)),arange()))) + (arange(1,numMUs)) - 0.4)],[time_spikes_AD(arange(dot((k - plotLengthWindows + 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp))),muActivation(arange(dot((k - plotLengthWindows + 1),ComParameters.nsamp) + 1,dot(k,ComParameters.nsamp)),arange()) / 7.5]]),h,ax,concat([time_DR(k - plotLengthWindows + 1),time_DR(k) + plotTimeRange / 2]),concat([0,1,1]))
            pause(ComParameters.nsamp / fs)
        tit=concat(['Block : ',num2str(k),' of ',num2str(nwin)])
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:182
        waitbar((k) / (nwin),f,tit)
        k=k + 1
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:184

    
    close_(f)
    for m in arange(1,size(signal.spikeTrains,2)).reshape(-1):
        signal.Dischargetimes[m]=find(signal.spikeTrains(arange(),m).T > 0)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:188
    
    __,dischargeRates=getFiringProperties(signal.spikeTrains.T,'outlierFlag',false,'flagPlotSpikeTrains',plotResults,'flagPlotDR',plotResults,'fs_MN',fs,'keepAllMUs',true,nargout=2)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:192
    # Store signal and activation variables
    signal.ct = copy(ct_dec)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:197
    signal.dischargeRates = copy(dischargeRates.T)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:198
    activation.muActivation = copy(muActivation)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:199
    activation.totalActivation = copy(totalActivation)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:200
    activation.ct = copy(ct_AD)
# ..\decompositionproject_python\decomposition_lib\simulate_online_decomposition.m:201
    return signal,decompParameters,activation
    
if __name__ == '__main__':
    pass
    