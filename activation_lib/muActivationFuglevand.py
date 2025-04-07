# Generated with SMOP  0.41-beta
from libsmop import *
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m

    
@function
def muActivationFuglevand(spikeTrains=None,Ap=None,tc=None,initialValues=None,varargin=None,*args,**kwargs):
    varargin = muActivationFuglevand.varargin
    nargin = muActivationFuglevand.nargin

    # muActivation calculates motor unit activation profiles.
# totalActivation corresponds to the sum (MU specific activation dynamics)
# Inputs:
# - spikeTrains: Matrix of spike trains (rows: time steps, columns: motor units)
# - Ap: Amplitude parameter
# - tc: Time constant parameter
# - initialValues: Neural activation history (if no history: zeros of size 2 x nMUs)
# Optional inputs (Name-Value pairs):
# - 'fs': Sampling frequency (default: 2048)
# - 'normalization': Type of normalization (default: 0)
# - 'twitchModel': Twitch model type (default: 'Fuglevand')
# - 'twitchSaturation': Enable twitch saturation (default: false)
# - 'satLevel': Saturation level (default: 0.3)
# - 'satSpeed': Saturation speed (default: 100)
    
    ## default values
    fs=2048
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:18
    normalization=0
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:19
    # todo improve saturation formulation
    numMUs=size(spikeTrains,2)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:21
    tetanusSaturation=copy(false)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:22
    satLevel=0.3
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:23
    satSpeed=100
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:24
    normSat=multiply(numMUs,satLevel)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:25
    ## Process optional input parameters using a loop and switch case
    for i in arange(1,length(varargin),2).reshape(-1):
        param=varargin[i]
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:29
        value=varargin[i + 1]
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:30
        if 'fs' == param:
            fs=copy(value)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:33
        else:
            if 'normalization' == param:
                normalization=copy(value)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:35
            else:
                if 'twitchSaturation' == param:
                    twitchSaturation=copy(value)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:37
                else:
                    if 'satLevel' == param:
                        satLevel=copy(value)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:39
                    else:
                        if 'satSpeed' == param:
                            satSpeed=copy(value)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:41
                        else:
                            error('Invalid optional input parameter: %s',param)
    
    ## initialize parameters
    muActivation=concat([[initialValues(concat([size(initialValues,1) - 1,size(initialValues,1)]),arange())],[zeros(size(spikeTrains))]])
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:48
    spikeTrains=concat([[zeros(2,size(spikeTrains,2))],[spikeTrains]])
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:50
    T=1 / fs
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:52
    ## normalization method
    if normalization == 1:
        updatedTwitchFactor=zeros(concat([1,size(spikeTrains,2)]))
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:56
        r=concat([0,0.4,0.75,1])
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:57
        f=concat([0,0.2,0.5,1])
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:58
        recCurve=fit(f.T,r.T,'poly2')
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:59
        normFactor=1.0 / numMUs
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:60
    else:
        if normalization == 2:
            # and Oye et al.
            m=0.045
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:63
            normFactor=1.0 / numMUs
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:64
            updatedTwitchFactor=zeros(concat([1,size(spikeTrains,2)]))
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:65
        else:
            if normalization == 3:
                normFactor=1.0 / numMUs
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:67
            else:
                if normalization == 4:
                    normFactor=1.0 / numMUs
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:69
                    Ap=atanh(Ap)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:70
    
    ## Start sample per sample decoding
    for n in arange(3,length(spikeTrains)).reshape(-1):
        #if sum(spikeTrains(n,:))>0 # to ensure that only after the first spike compute the models
        if normalization == 0:
            muActivation[n,arange()]=multiply((dot(2,exp(- T / tc))).T,muActivation(n - 1,arange())) - multiply((exp((dot(- 2,T)) / tc)).T,muActivation(n - 2,arange())) + multiply((multiply(((dot(Ap,T)) / tc),exp(1 - T / (tc)))).T,spikeTrains(n - 1,arange()))
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:77
        else:
            if normalization == 1:
                # updateTwitchFactor updates the twitch only when a new twitch appears
                updatedTwitchFactor[logical(spikeTrains(n - 1,arange()))]=dot(1.5,recCurve(abs(force(n)))) + 0.5
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:81
                muActivation[n,arange()]=multiply((dot(2,exp(- T / tc))).T,muActivation(n - 1,arange())) - multiply((exp((dot(- 2,T)) / tc)).T,muActivation(n - 2,arange())) + multiply(multiply((multiply(((dot(dot((normFactor),Ap),T)) / tc),exp(1 - T / (tc)))).T,updatedTwitchFactor),spikeTrains(n - 1,arange()))
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:82
                #                     neuralActivation(n,:) = (2*exp(-T./tc))'.*neuralActivation(n-1,:)  -  (exp((-2*T)./tc))'.*neuralActivation(n-2,:)  ...
            #                         + ((((normFactor(i)*recCurve(force(n)))*Ap*T)./tc).*exp(1-T./(tc)))'.*excitation(n-1,:);
            else:
                if normalization == 2:
                    updatedTwitchFactor[logical(spikeTrains(n - 1,arange()))]=((log(dot(100,abs(force(n - 1)))) / m) - 2.1178) / 100
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:88
                    if updatedTwitchFactor(logical(spikeTrains(n - 1,arange()))) >= 1:
                        updatedTwitchFactor[logical(spikeTrains(n - 1,arange()))]=1
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:90
                    muActivation[n,arange()]=multiply((dot(2,exp(- T / tc))).T,muActivation(n - 1,arange())) - multiply((exp((dot(- 2,T)) / tc)).T,muActivation(n - 2,arange())) + multiply(multiply((multiply(((dot(dot((normFactor),Ap),T)) / tc),exp(1 - T / (tc)))).T,updatedTwitchFactor),spikeTrains(n - 1,arange()))
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:92
                else:
                    if normalization >= 3:
                        muActivation[n,arange()]=multiply((dot(2,exp(- T / tc))).T,muActivation(n - 1,arange())) - multiply((exp((dot(- 2,T)) / tc)).T,muActivation(n - 2,arange())) + multiply((multiply(((dot(dot(dot(((normFactor)),abs(force(n) - force(n - 1))),Ap),T)) / tc),exp(1 - T / (tc)))).T,spikeTrains(n - 1,arange()))
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:95
        if tetanusSaturation:
            neuralActivationTransf[n,arange()]=dot(normSat,(1 - dot(1,exp(dot(- satSpeed,(neuralActivation(n,arange()))))))) / (1 + dot(1,exp(dot(- satSpeed,(neuralActivation(n,arange()))))))
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:99
        #end
    
    if tetanusSaturation:
        muActivation=copy(neuralActivationTransf)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:104
    
    muActivation[arange(1,2),arange()]=[]
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:106
    
    totalActivation=sum(muActivation,2)
# ../decompositionproject_python/activation_lib\muActivationFuglevand.m:107
    return muActivation,totalActivation
    
if __name__ == '__main__':
    pass
    