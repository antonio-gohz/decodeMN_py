# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\remoutliers.m

    ###########################################################################
# To remove high values of DR that we consider as outliers until we reach a
# low coefficient of variation of the discharge rates
    
    # Input: 
#   PulseTold = previous pulse trains of the motor units
#   Distimeold = previous discharge times of the motor units
#   thresh = Threshold for the coefficent of variation of the discharge
#   rate that we want to reach
#   fsamp = sampling frequency
    
    # Output:
#   Distime = new discharge times of the motor units
    
    ###########################################################################
    
    
@function
def remoutliers(pulseT=None,distime=None,thresh=None,fsamp=None,*args,**kwargs):
    varargin = remoutliers.varargin
    nargin = remoutliers.nargin

    for nMU in arange(1,length(distime)).reshape(-1):
        DR=1.0 / (diff(distime[nMU]) / fsamp)
# ..\decompositionproject_python\decomposition_lib\remoutliers.m:22
        k=1
# ..\decompositionproject_python\decomposition_lib\remoutliers.m:24
        while (std(DR) / mean(DR)) > thresh and k < 30:

            k=k + 1
# ..\decompositionproject_python\decomposition_lib\remoutliers.m:26
            thres=mean(DR) + dot(3,std(DR))
# ..\decompositionproject_python\decomposition_lib\remoutliers.m:27
            idx=find(DR > thres)
# ..\decompositionproject_python\decomposition_lib\remoutliers.m:28
            if logical_not(isempty(idx)):
                for i in arange(1,length(idx)).reshape(-1):
                    if pulseT(nMU,distime[nMU](idx(i))) < pulseT(nMU,distime[nMU]((idx(i) + 1))):
                        idxdel[i]=idx(i)
# ..\decompositionproject_python\decomposition_lib\remoutliers.m:32
                    else:
                        idxdel[i]=idx(i) + 1
# ..\decompositionproject_python\decomposition_lib\remoutliers.m:34
                distime[nMU][idxdel]=[]
# ..\decompositionproject_python\decomposition_lib\remoutliers.m:37
                clearvars('idxdel')
            DR=1.0 / (diff(distime[nMU]) / fsamp)
# ..\decompositionproject_python\decomposition_lib\remoutliers.m:40

    