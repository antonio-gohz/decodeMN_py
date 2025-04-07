# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m

    
@function
def getonlineparameters(eSIG=None,MUFilters=None,fsamp=None,*args,**kwargs):
    varargin = getonlineparameters.varargin
    nargin = getonlineparameters.nargin

    ReSIG=dot(eSIG,eSIG.T) / length(eSIG)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:3
    iReSIGt=pinv(ReSIG)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:4
    nMU=1
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:6
    for i in arange(1,size(MUFilters,2)).reshape(-1):
        PulseT=dot((dot(MUFilters(arange(),nMU).T,iReSIGt)),eSIG)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:8
        PulseT=multiply(PulseT,abs(PulseT))
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:10
        __,spikes=findpeaks(PulseT,'MinPeakDistance',round(dot(fsamp,0.02)),nargout=2)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:11
        if length(spikes) > 2:
            L,C=kmeans(PulseT(spikes).T,2,nargout=2)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:13
            __,idx=max(C,nargout=2)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:14
            peakspikes=spikes(L == idx)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:15
            peaksnoise=setdiff(spikes,peakspikes)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:16
            peakspikes[PulseT(peakspikes) > mean(PulseT(peakspikes)) + dot(3,std(PulseT(peakspikes)))]=[]
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:17
            norm[i]=mean(maxk(PulseT(peakspikes),10))
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:18
            PulseT=PulseT / mean(maxk(PulseT(spikes),10))
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:19
            __,centroid(i,1)=kmeans(PulseT(peaksnoise).T,1,nargout=2)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:20
            __,centroid(i,2)=kmeans(PulseT(peakspikes).T,1,nargout=2)
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:21
            nMU=nMU + 1
# ..\decompositionproject_python\decomposition_lib\getonlineparameters.m:22
    