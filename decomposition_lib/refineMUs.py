# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\refineMUs.m

    ###########################################################################
# To reapply the filters of each motor unit over the EMG entire signal
    
    # Input: 
#   EMG = raw EMG signal
#   EMGmask = flagged EMG channels with artifacts or low signal to noise
#   ratio
#   PulseTold = previous pulse trains of the motor units
#   Distimeold = previous discharge times of the motor units
#   fsamp = sampling frequency
    
    # Output:
#   PulseT = new pulse trains of the motor units
#   Distime = new discharge times of the motor units
    
    ###########################################################################
    
    
@function
def refineMUs(EMG=None,EMGmask=None,PulseTold=None,Distimeold=None,fsamp=None,*args,**kwargs):
    varargin = refineMUs.varargin
    nargin = refineMUs.nargin

    f=waitbar(0,'Refining MU pulse trains - Signal extension')
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:20
    EMG[EMGmask == 1,arange()]=[]
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:22
    nbextchan=1500
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:23
    exFactor=round(nbextchan / size(EMG,1))
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:24
    eSIG=extend(EMG,exFactor)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:25
    ReSIG=dot(eSIG,eSIG.T) / length(eSIG)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:26
    iReSIGt=pinv(ReSIG)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:27
    PulseT=zeros(size(PulseTold))
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:29
    Distime=cell(1,size(PulseTold,1))
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:30
    x=1 / size(PulseTold,1)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:32
    # Recalculate MUfilters
    for i in arange(1,size(PulseTold,1)).reshape(-1):
        Distimeold[i][PulseTold(i,Distimeold[i]) > mean(PulseTold(i,Distimeold[i])) + dot(3,std(PulseTold(i,Distimeold[i])))]=[]
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:36
        MUFilters=sum(eSIG(arange(),Distimeold[i]),2)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:37
        IPTtmp=dot((dot(MUFilters.T,iReSIGt)),eSIG)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:38
        PulseT[i,arange()]=IPTtmp(arange(1,size(EMG,2)))
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:39
        PulseT[i,arange()]=multiply(abs(PulseT(i,arange())),PulseT(i,arange()))
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:40
        __,spikes=findpeaks(PulseT(i,arange()),'MinPeakDistance',round(dot(fsamp,0.02)),nargout=2)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:41
        PulseT[i,arange()]=PulseT(i,arange()) / mean(maxk(PulseT(i,spikes),10))
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:42
        L,C=kmeans(PulseT(i,spikes).T,2,nargout=2)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:43
        __,idx=max(C,nargout=2)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:44
        Distime[i]=spikes(L == idx)
# ..\decompositionproject_python\decomposition_lib\refineMUs.m:45
        waitbar(dot(x,i),f,concat(['Refining MU#',num2str(i),' out of ',num2str(size(PulseTold,1)),'MUs']))
    
    close_(f)