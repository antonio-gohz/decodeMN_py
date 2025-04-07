# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m

    
@function
def reuseMUFilters(EMGs=None,decompParameters=None,varargin=None,*args,**kwargs):
    varargin = reuseMUFilters.varargin
    nargin = reuseMUFilters.nargin

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
    if size(EMGs,1) > size(EMGs,2):
        EMGs=EMGs.T
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:17
    
    # default parameters
    fs=2048
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:22
    for i in arange(1,length(varargin),2).reshape(-1):
        param=varargin[i]
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:25
        value=varargin[i + 1]
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:26
        if 'fs' == param:
            fs=copy(value)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:29
        else:
            error('Invalid optional input parameter: %s',param)
    
    tanh_denoise=3
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:34
    # whitened MU filter * whitened data
    tic
    X=extend(EMGs(decompParameters.EMGmask,arange()),decompParameters.extensionfactor)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:37
    #X = decompParameters.whitMat*eSIG;
#normtanh =  tanh_denoise*std(X,[],2);
#if ~isempty(decompParameters.normtanh)
    #X = decompParameters.normtanh.*tanh(X./decompParameters.normtanh);
#    X = normtanh.*tanh(X./normtanh);
    
    #end
    
    # PulseT = (decompParameters.MUFilters' * X).*abs(decompParameters.MUFilters' * X); # 4a: Estimate the i source
    PulseT=(dot(dot(decompParameters.MUFilters.T,decompParameters.iReSIGt),X))
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:47
    PulseT=multiply(PulseT,abs(PulseT))
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:48
    PulseT=PulseT / decompParameters.normIPT.T
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:49
    PulseT=tanh(PulseT)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:50
    spikes1,__=islocalmax(PulseT,2,'MinSeparation',round(dot(fs,0.02)),nargout=2)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:52
    spikeTrains=multiply(PulseT,spikes1) >= abs(diff(decompParameters.centroid,1,2)) / 2 + min(decompParameters.centroid,[],2)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:54
    ct=copy(toc)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:55
    __,__,Distime=getFiringProperties(spikeTrains.T,'flagPlotDR',0,'flagPlotSpikeTrains',0,nargout=3)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:56
    SIL=getSIL_PNR(spikeTrains,PulseT,0,decompParameters.centroid)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:59
    PulseT,Distime,idsnewVec=remduplicates(PulseT,Distime,Distime,round(fs / 40),0.00025,0.3,fs,SIL,nargout=3)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:60
    decompParameters.MUFilters = copy(decompParameters.MUFilters(arange(),idsnewVec))
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:62
    decompParameters.normIPT = copy(decompParameters.normIPT(idsnewVec))
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:63
    decompParameters.centroid = copy(decompParameters.centroid(idsnewVec,arange()))
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:64
    #decompParameters.normtanh =normtanh;
    
    signal.spikeTrains = copy(spikeTrains(idsnewVec,arange()))
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:67
    signal.Pulsetrain = copy(PulseT)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:68
    signal.SIL = copy(SIL(idsnewVec))
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:69
    signal.Dischargetimes = copy(Distime)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:70
    signal.idsnew = copy(idsnewVec)
# ..\decompositionproject_python\decomposition_lib\reuseMUFilters.m:71
    
    #spikeTrains = (abs(PulseT .* spikes1 - min(decompParameters.centroid,[],2)) > abs(PulseT .* spikes1 - max(decompParameters.centroid,[],2)));
    
    #Dischargetimes
    
    # dewhitened filter * iRSe * extended data
# MUFiltersdeWhit = dewhiteningMatrix * MUFilters;
# 
# ReSIG = eSIG*eSIG'/length(eSIG);
# iReSIGt = pinv(ReSIG);
# icasig2 = (MUFiltersdeWhit' * iReSIGt*eSIG).*abs(MUFiltersdeWhit' * iReSIGt*eSIG); # 4a: Estimate the i source
#  hold on; plot(icasig2(1,:))
# 
# PulseTtmp = (MUFilters(:,nMU)'*iReSIGt)*eSIG;
# PulseT(nMU,:) = PulseTtmp(1:size(EMG,2));
# PulseT(nMU,:) = PulseT(nMU,:) .* abs(PulseT(nMU,:));
# icasig = (MUFiltersdeWhit' * X).*abs(MUFiltersdeWhit' * X); # 4a: Estimate the i source
    