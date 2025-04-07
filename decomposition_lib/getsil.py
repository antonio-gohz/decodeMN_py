# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\getsil.m

    
@function
def getsil(PulseT=None,fsamp=None,*args,**kwargs):
    varargin = getsil.varargin
    nargin = getsil.nargin

    __,spikes=findpeaks(PulseT,'MinPeakDistance',round(dot(fsamp,0.02)),nargout=2)
# ..\decompositionproject_python\decomposition_lib\getsil.m:3
    
    PulseT=PulseT / mean(maxk(PulseT(spikes),10))
# ..\decompositionproject_python\decomposition_lib\getsil.m:4
    
    L,C,sumd,D=kmeans(PulseT(spikes).T,2,nargout=4)
# ..\decompositionproject_python\decomposition_lib\getsil.m:5
    
    __,idx2=max(C,nargout=2)
# ..\decompositionproject_python\decomposition_lib\getsil.m:6
    
    within=sumd(idx2)
# ..\decompositionproject_python\decomposition_lib\getsil.m:7
    between=sum(D(L == idx2,setdiff(concat([1,2]),idx2)))
# ..\decompositionproject_python\decomposition_lib\getsil.m:8
    sil=(between - within) / max(concat([within,between]))
# ..\decompositionproject_python\decomposition_lib\getsil.m:9
    