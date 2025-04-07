# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m

    
@function
def maximizeSIL_online(X=None,SIL=None,w=None,fsamp=None,showPlots=None,h=None,ax=None,normIPT=None,centroids=None,*args,**kwargs):
    varargin = maximizeSIL_online.varargin
    nargin = maximizeSIL_online.nargin

    ###########################################################################
# Optimization loop of the MU filter to minimize the Silhoutte
    
    # Input: 
#   w = initial weigths
#   X = whitened signal
#   SIL = Silhoutte
#   fsamp = sampling frequency
    
    # Output:
#   wlast = new weigths (MU filter)
#   spikeslast = discharge times of the motor unit
#   SILlast = coefficient of varation of the inter spike intervals
    
    ###########################################################################
    k=1
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:19
    SILlast=copy(eps)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:20
    xX=linspace(0,length(X) / fsamp,length(X))
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:21
    xW=arange(1,length(w))
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:22
    spikes=zeros(1,size(X,2))
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:23
    ipt=zeros(size(X,2),1)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:24
    while SIL > SILlast:

        SILlast=copy(SIL)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:28
        wlast=copy(w)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:29
        normIPTlast=copy(normIPT)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:30
        centroidslast=copy(centroids)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:31
        icasig=(multiply((dot(w.T,X)),abs(dot(w.T,X))))
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:34
        ipt=icasig / normIPT.T
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:35
        ipt=tanh(ipt)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:36
        spikes,__=islocalmax(ipt.T,1,'MinSeparation',round(dot(fsamp,0.02)),nargout=2)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:37
        spikes=multiply(ipt.T,spikes) >= (abs(diff(centroids,1,2)) / 2 + min(centroids,[],2)).T
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:38
        # at least 5 pps from the first spike to the last detected
        sufficientSpikes=find(spikes)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:42
        if logical_not(isempty(sufficientSpikes)):
            sufficientSpikes=round(dot(5.0,diff(concat([sufficientSpikes(1),sufficientSpikes(end())]))) / fsamp)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:44
            if sum(spikes) > sufficientSpikes:
                SIL,__,centroids=getsil(ipt,nargout=3)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:46
                normIPT=mean(maxk(icasig(spikes),10))
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:48
                # better to calculate SIL without tanh so the original is lower an
            # it can be improved
                ipt=icasig / normIPT.T
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:52
                ipt=tanh(ipt)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:53
                spikes=multiply(ipt.T,spikes) >= (abs(diff(centroids,1,2)) / 2 + min(centroids,[],2)).T
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:54
                SIL,__,centroids=getSIL_PNR(spikes,ipt,nargout=3)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:55
                w=sum(X(arange(),spikes),2)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:59
                w=w / norm(w)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:60
                if showPlots:
                    h,ax=livePlots.updatePlot(concat([[cellarray([xX,ipt / max(ipt)])],[cellarray([xX(spikes.T),ipt(spikes.T) / max(ipt)])],[cellarray([xW,w.T])]]),h,ax,nargout=2)
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:62
                    title(ax(1),'2 Refinement SIL: ' + SILlast + ' It: ' + k)
        k=k + 1
# ..\decompositionproject_python\decomposition_lib\maximizeSIL_online.m:68

    