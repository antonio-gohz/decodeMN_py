# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\realignspikes.m

    ###########################################################################
# Function to realign each discharge time with the maxima of each MU action
# potential (identified over the channel with the highest peak-to-peak
# amplitude)
    
    # Input: 
#   EMG = raw EMG signal
#   EMGmask = flagged EMG channels with artifacts or low signal to noise
#   ratio
#   coordinates = x and y coordinates of each EMG channel over the grid of
#   electrodes
#   Distimeold = previous discharge times for spike trigger identification
#   of MUAP
#   fsamp = sampling frequency
#   win = duration of the window for spike trigger averaging
    
    # Output:
#   Distime = new discharge times of the motor units
    
    ###########################################################################
    
    
@function
def realignspikes(EMG=None,EMGmask=None,coordinates=None,Distimeold=None,fsamp=None,win=None,*args,**kwargs):
    varargin = realignspikes.varargin
    nargin = realignspikes.nargin

    EMG[EMGmask == 1,arange()]=[]
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:24
    EMG2=cell(max(coordinates(arange(),1)),max(coordinates(arange(),2)))
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:25
    # reorganize the EMG data in 2d cells before performing a double differential
# montage for each channel
    for i in arange(1,size(EMG,1)).reshape(-1):
        EMG2[coordinates(i,1),coordinates(i,2)]=EMG(i,arange())
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:29
    
    # To generate the double differential EMG signals
    ch=1
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:33
    for c in arange(1,size(EMG2,2)).reshape(-1):
        for r in arange(1,size(EMG2,1) - 2).reshape(-1):
            if logical_not(isempty(EMG2[r,c])) and logical_not(isempty(EMG2[r + 1,c])) and logical_not(isempty(EMG2[r + 2,c])):
                EMGdiff[ch,arange()]=(EMG2[r,c] - EMG2[r + 1,c]) - (EMG2[r + 1,c] - EMG2[r + 2,c])
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:37
                ch=ch + 1
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:38
    
    Distime=cell(1,length(Distimeold))
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:43
    for i in arange(1,length(Distimeold)).reshape(-1):
        p2p=zeros(1,size(EMGdiff,1))
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:45
        for l in arange(1,size(EMGdiff,1)).reshape(-1):
            temp=cutMUAP(Distimeold[i],round(dot(win,fsamp)),EMGdiff(l,arange()))
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:47
            p2p[l]=max(mean(temp,1)) - min(mean(temp,1))
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:48
        __,idxEMG=max(p2p,nargout=2)
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:50
        Distime[i]=alignMUAP(Distimeold[i],round(dot(win,fsamp)),EMGdiff(idxEMG,arange()))
# ..\decompositionproject_python\decomposition_lib\realignspikes.m:51
        clearvars('temp')
    
    return Distime
    
if __name__ == '__main__':
    pass
    