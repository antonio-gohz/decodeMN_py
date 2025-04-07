# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m

    
@function
def getSIL_PNR(spikeTrains=None,IPTs=None,plotFlag=None,centroid=None,*args,**kwargs):
    varargin = getSIL_PNR.varargin
    nargin = getSIL_PNR.nargin

    #UNTITLED Summary of this function goes here
#   Detailed explanation goes here
    
    centroidsFlag=copy(false)
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:5
    if nargin < 3:
        plotFlag=0
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:8
    
    if nargin < 4:
        centroidsFlag=copy(true)
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:11
        centroids=zeros(concat([1,2]))
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:12
    
    spikeTrains=logical(spikeTrains)
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:15
    if size(spikeTrains,1) > size(spikeTrains,2):
        spikeTrains=spikeTrains.T
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:18
    
    if size(IPTs,1) > size(IPTs,2):
        IPTs=IPTs.T
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:22
    
    sil=zeros(1,size(spikeTrains,1))
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:24
    pnr=zeros(1,size(spikeTrains,1))
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:25
    for i in arange(1,size(spikeTrains,1)).reshape(-1):
        if centroidsFlag:
            centroids[1]=mean(IPTs(i,logical_not(spikeTrains(i,arange()))))
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:28
            centroids[2]=mean(IPTs(i,spikeTrains(i,arange())))
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:29
        else:
            centroids[1]=min(centroid(i,arange()))
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:31
            centroids[2]=max(centroid(i,arange()))
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:32
        # euclidean distance within spike cluster
        within=dot((IPTs(i,spikeTrains(i,arange())) - centroids(2)),(IPTs(i,spikeTrains(i,arange())) - centroids(2)).T)
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:35
        between=dot((IPTs(i,spikeTrains(i,arange())) - centroids(1)),(IPTs(i,spikeTrains(i,arange())) - centroids(1)).T)
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:37
        sil[i]=(between - within) / max(concat([within,between]))
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:38
        pnr[i]=dot(10,log10(mean(IPTs(i,spikeTrains(i,arange())) ** 2) / mean(IPTs(i,logical_not(spikeTrains(i,arange()))) ** 2)))
# ..\decompositionproject_python\decomposition_lib\getSIL_PNR.m:39
        if plotFlag:
            figure
            plot(IPTs(i,arange()))
            hold('on')
            plot(find(spikeTrains(i,arange())),IPTs(i,spikeTrains(i,arange())),'or')
            #silhouette(IPTs(i,:)',spikeTrains(i,:))
            title('SIL: ' + sil(i) + ', PNR: ' + pnr(i) + ' dB')
    
    return sil,pnr,centroids
    
if __name__ == '__main__':
    pass
    