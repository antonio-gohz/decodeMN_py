# Generated with SMOP  0.41-beta
from libsmop import *
# ../decompositionproject_python/activation_lib\activationDynamics.m

    
@function
def activationDynamics(excitation=None,dischargeRate=None,t=None,Ap=None,tc=None,thr=None,twitchModel=None,satLevel=None,tetParam=None,tetanusModel=None,c1=None,c2=None,shapeFactor_=None,activationScale_=None,force=None,muscle=None,subject=None,position=None,activation=None,normalization=None,*args,**kwargs):
    varargin = activationDynamics.varargin
    nargin = activationDynamics.nargin

    neuralActivation=zeros(size(excitation))
# ../decompositionproject_python/activation_lib\activationDynamics.m:4
    neuralMatrix=zeros(max(sum(excitation)),size(excitation,2))
# ../decompositionproject_python/activation_lib\activationDynamics.m:5
    n_spike=0
# ../decompositionproject_python/activation_lib\activationDynamics.m:6
    firstSpike=0
# ../decompositionproject_python/activation_lib\activationDynamics.m:7
    t0s=NaN(max(sum(excitation)),size(excitation,2))
# ../decompositionproject_python/activation_lib\activationDynamics.m:8
    
    
    T=t(2) - t(1)
# ../decompositionproject_python/activation_lib\activationDynamics.m:11
    normSat=multiply(1 / size(excitation,2),satLevel)
# ../decompositionproject_python/activation_lib\activationDynamics.m:12
    
    ## Tetanus parameters
    if tetanusModel:
        #tetParam=100; # the bigger ther smoother saturation
        normSat=multiply(1 / size(excitation,2),satLevel)
# ../decompositionproject_python/activation_lib\activationDynamics.m:18
    
    ## Raikova parameters
    if strcmp(twitchModel,'Raikova'):
        k=log(2) / (thr - tc - multiply(tc,log(thr / tc)))
# ../decompositionproject_python/activation_lib\activationDynamics.m:22
        m=multiply(k,tc)
# ../decompositionproject_python/activation_lib\activationDynamics.m:23
        p=multiply(Ap,exp(multiply(multiply(- k,tc),(log(tc) - 1))))
# ../decompositionproject_python/activation_lib\activationDynamics.m:24
    
    ## CEINMS variables
    beta1_=c1 + c2
# ../decompositionproject_python/activation_lib\activationDynamics.m:27
    beta2_=multiply(c1,c2)
# ../decompositionproject_python/activation_lib\activationDynamics.m:28
    alpha_=1 + beta1_ + beta2_
# ../decompositionproject_python/activation_lib\activationDynamics.m:29
    neuralActivation_=copy(neuralActivation)
# ../decompositionproject_python/activation_lib\activationDynamics.m:30
    
    tStamps=extractTimeStamps(muscle,subject,position,activation)
# ../decompositionproject_python/activation_lib\activationDynamics.m:33
    numMUs=[]
# ../decompositionproject_python/activation_lib\activationDynamics.m:34
    fs=2048
# ../decompositionproject_python/activation_lib\activationDynamics.m:35
    try:
        for i in arange(1,size(tStamps,1)).reshape(-1):
            numMUs=concat([numMUs,sum(sum(excitation(arange(floor(dot(fs,tStamps(i,1))),floor(dot(fs,tStamps(i,2)))),arange())) > 0)])
# ../decompositionproject_python/activation_lib\activationDynamics.m:38
    finally:
        pass
    
    # This is a small band aid very large differences of the number of MUs
    #make the normalization amplityde very different e.g. 1/5 vs 1/30 MUs
    # this happens a lot in the soleau theregore limit the numMUs to a
    # minimum of 10 MUs
    numMUs[numMUs < 10]=10
# ../decompositionproject_python/activation_lib\activationDynamics.m:47
    numMUs[numMUs > 30]=30
# ../decompositionproject_python/activation_lib\activationDynamics.m:48
    if normalization == 1:
        updatedTwitchFactor=zeros(concat([1,size(excitation,2)]))
# ../decompositionproject_python/activation_lib\activationDynamics.m:52
        #                 figure;
        r=concat([0,0.4,0.75,1])
# ../decompositionproject_python/activation_lib\activationDynamics.m:54
        f=concat([0,0.2,0.5,1])
# ../decompositionproject_python/activation_lib\activationDynamics.m:55
        recCurve=fit(f.T,r.T,'poly2')
# ../decompositionproject_python/activation_lib\activationDynamics.m:56
        # compare with
#Determine percentage of activation (duchateau and Enoka 2021) TA
#         f=0:0.01:1;
#             hold on,plot(f,recCurve(f));
# m = 0.0387;
# proportion = ((log(100*f)/m) - 2.1178)./100;
# 
# hold on, plot (f,proportion)
# 
# m = 0.0446;
# proportion = ((log(100*f)/m) - 2.1178)./100;
# hold on, plot (f,proportion)
# 
# rec_sol_MUs = [8,5,3,4,2,2,5,6,4,3];
# rec_sol_cumMUs = cumsum(rec_sol_MUs)/sum(rec_sol_MUs);
# rec_sol_activation = 0.025:.10:1;
# hold on, plot(rec_sol_cumMUs,rec_sol_activation)
# xlabel('Torque # MVC')
# ylabel('MU rec(0-1)')
# 
# legend('Old generic Curve', 'TA curve', 'Soleus', 'Soleus, (literature)') 
#     switch muscle
#         case 'TA',    totalMUs = 350;
#         case 'PERTERT',  totalMUs = 200; #check
#         case 'GASlat',  totalMUs = 500;
#         case 'GASmed',  totalMUs = 400;
#         case 'PERLONGUS',  totalMUs = 300; #check
#         case 'SOL',  totalMUs = 900;
#     end
        normFactor=1.0 / numMUs
# ../decompositionproject_python/activation_lib\activationDynamics.m:87
    else:
        if normalization == 2:
            # and Oye et al.
        #numMUs=ones(size(numMUs)); # test get rid of it
            m=0.045
# ../decompositionproject_python/activation_lib\activationDynamics.m:92
            # WRONG : I used a scaled version of rafa for m //
#         switch muscle   obtained with script finalPlots, based on TA and SOL refereces below: 
#             case 'TA',    m = 0.0387;  totalMUs = 350; # Enoka to reach 1 at 50# recruited
#             case 'PERTERT', m = 0.0386; totalMUs = 200;
#             case 'GASlat',  m = 0.0366; totalMUs = 500;
#             case 'GASmed',  m = 0.0371; totalMUs = 400;
#             case 'PERLONGUS', m = 0.0396; totalMUs = 300;
#             case 'SOL',  m = 0.0446; totalMUs = 900; # oye et al to reach 1 of act at 95# recruited
#         end
            normFactor=1.0 / numMUs
# ../decompositionproject_python/activation_lib\activationDynamics.m:102
            updatedTwitchFactor=zeros(concat([1,size(excitation,2)]))
# ../decompositionproject_python/activation_lib\activationDynamics.m:104
        else:
            if normalization == 3:
                normFactor=1.0 / numMUs
# ../decompositionproject_python/activation_lib\activationDynamics.m:108
            else:
                if normalization == 4:
                    normFactor=1.0 / numMUs
# ../decompositionproject_python/activation_lib\activationDynamics.m:110
                    Ap=atanh(Ap)
# ../decompositionproject_python/activation_lib\activationDynamics.m:111
    
    i=1
# ../decompositionproject_python/activation_lib\activationDynamics.m:120
    
    for n in arange(3,length(excitation)).reshape(-1):
        if sum(excitation(n,arange())) > 0:
            n_spike=n_spike + 1
# ../decompositionproject_python/activation_lib\activationDynamics.m:125
            t0=t(1) + dot((n - 1),T)
# ../decompositionproject_python/activation_lib\activationDynamics.m:126
            firstSpike=1
# ../decompositionproject_python/activation_lib\activationDynamics.m:127
            t0s[n_spike,excitation(n,arange()) > 0]=t0
# ../decompositionproject_python/activation_lib\activationDynamics.m:128
        # i is the counter of the ramps as each one has a different
        # decomposition, i.e., a different number of MUs and hence
        # different normFactor
        # if a new ramp (n>floor(tStamps(i+1,1))) and i is not bigger than
        # the total amount of ramps (size(tStamps,1))
        if i < size(tStamps,1) and n > floor(dot(fs,tStamps(i + 1,1))):
            i=i + 1
# ../decompositionproject_python/activation_lib\activationDynamics.m:137
        #         if n_spike >= 554  # just for debugging 
#             neuralMatrix(553,23)'
#         end
        if firstSpike:
            if strcmp(twitchModel,'Raikova'):
                neuralMatrix[arange(1,n_spike),arange()]=(multiply(multiply(p,round(t(n) - t0s(arange(1,n_spike),arange()),8).T ** m),exp(multiply(- k,round(t(n) - t0s(arange(1,n_spike),arange()),8).T)))).T
# ../decompositionproject_python/activation_lib\activationDynamics.m:145
                neuralActivation[n,arange()]=sum(neuralMatrix(arange(1,n_spike),arange()),1,'omitnan')
# ../decompositionproject_python/activation_lib\activationDynamics.m:146
            else:
                if strcmp(twitchModel,'Fuglevand'):
                    if normalization == 0:
                        neuralActivation[n,arange()]=multiply((dot(2,exp(- T / tc))).T,neuralActivation(n - 1,arange())) - multiply((exp((dot(- 2,T)) / tc)).T,neuralActivation(n - 2,arange())) + multiply((multiply(((dot(Ap,T)) / tc),exp(1 - T / (tc)))).T,excitation(n - 1,arange()))
# ../decompositionproject_python/activation_lib\activationDynamics.m:149
                    else:
                        if normalization == 1:
                            # updateTwitchFactor updates the twitch only when a new
                    # a new twitch appears
                            updatedTwitchFactor[logical(excitation(n - 1,arange()))]=dot(1.5,recCurve(abs(force(n)))) + 0.5
# ../decompositionproject_python/activation_lib\activationDynamics.m:154
                            neuralActivation[n,arange()]=multiply((dot(2,exp(- T / tc))).T,neuralActivation(n - 1,arange())) - multiply((exp((dot(- 2,T)) / tc)).T,neuralActivation(n - 2,arange())) + multiply(multiply((multiply(((dot(dot((normFactor(i)),Ap),T)) / tc),exp(1 - T / (tc)))).T,updatedTwitchFactor),excitation(n - 1,arange()))
# ../decompositionproject_python/activation_lib\activationDynamics.m:155
                            # practically the same the one above waits for next twitch 
                    #                     neuralActivation(n,:) = (2*exp(-T./tc))'.*neuralActivation(n-1,:)  -  (exp((-2*T)./tc))'.*neuralActivation(n-2,:)  ...
#                         + ((((normFactor(i)*recCurve(force(n)))*Ap*T)./tc).*exp(1-T./(tc)))'.*excitation(n-1,:);
                        else:
                            if normalization == 2:
                                #normalizing according to recruitmentc coding
                                updatedTwitchFactor[logical(excitation(n - 1,arange()))]=((log(dot(100,abs(force(n - 1)))) / m) - 2.1178) / 100
# ../decompositionproject_python/activation_lib\activationDynamics.m:162
                                if updatedTwitchFactor(logical(excitation(n - 1,arange()))) >= 1:
                                    updatedTwitchFactor[logical(excitation(n - 1,arange()))]=1
# ../decompositionproject_python/activation_lib\activationDynamics.m:164
                                #normalizing according to rate coding this might be
                    #wrong 
#                     if ~isnan(dischargeRate(n-1))
#                     updatedTwitchFactor(logical(excitation(n-1,:)))= dischargeRate(n-1).*updatedTwitchFactor(logical(excitation(n-1,:)));
#                     end
# if sum(~isnan(dischargeRate(n-1,find(logical(excitation(n-1,:))))))>=1
#     updatedTwitchFactor(logical(excitation(n-1,:)))= dischargeRate(n-1,find(logical(excitation(n-1,:)))).*...
#         updatedTwitchFactor(logical(excitation(n-1,:)));
# end
                                neuralActivation[n,arange()]=multiply((dot(2,exp(- T / tc))).T,neuralActivation(n - 1,arange())) - multiply((exp((dot(- 2,T)) / tc)).T,neuralActivation(n - 2,arange())) + multiply(multiply((multiply(((dot(dot((normFactor(i)),Ap),T)) / tc),exp(1 - T / (tc)))).T,updatedTwitchFactor),excitation(n - 1,arange()))
# ../decompositionproject_python/activation_lib\activationDynamics.m:177
                            else:
                                if normalization >= 3:
                                    neuralActivation[n,arange()]=multiply((dot(2,exp(- T / tc))).T,neuralActivation(n - 1,arange())) - multiply((exp((dot(- 2,T)) / tc)).T,neuralActivation(n - 2,arange())) + multiply((multiply(((dot(dot(dot(((normFactor(i))),abs(force(n) - force(n - 1))),Ap),T)) / tc),exp(1 - T / (tc)))).T,excitation(n - 1,arange()))
# ../decompositionproject_python/activation_lib\activationDynamics.m:180
                else:
                    if strcmp(twitchModel,'CEINMS'):
                        neuralActivation_[n,arange()]=dot(alpha_,excitation(n,arange())) - dot(beta1_,neuralActivation_(n - 1,arange())) - dot(beta2_,neuralActivation_(n - 2,arange()))
# ../decompositionproject_python/activation_lib\activationDynamics.m:184
                        neuralActivation_[n,arange()]=dot(activationScale_,(exp(dot(shapeFactor_,neuralActivation_(n,arange()))) - 1)) / (exp(shapeFactor_) - 1)
# ../decompositionproject_python/activation_lib\activationDynamics.m:186
                    else:
                        if strcmp(twitchModel,'Piecewise'):
                            rLessTc,cLessTc=find(t(n) - t0s(arange(1,n_spike),arange()) < tc,nargout=2)
# ../decompositionproject_python/activation_lib\activationDynamics.m:188
                            rInTcThr,cInTcThr=find((multiply((t(n) - t0s(arange(1,n_spike),arange()) > tc),(t(n) - t0s(arange(1,n_spike),arange()) < dot(2,(thr - tc)) + tc))) > 0,nargout=2)
# ../decompositionproject_python/activation_lib\activationDynamics.m:189
                            if n_spike >= 2:
                                rMoreThr,cMoreThr=find(t(n) - t0s(arange(1,n_spike),arange()) > dot(2,(thr - tc)) + tc,nargout=2)
# ../decompositionproject_python/activation_lib\activationDynamics.m:191
                                neuralMatrix[rMoreThr,cMoreThr]=0
# ../decompositionproject_python/activation_lib\activationDynamics.m:192
                            neuralMatrix[rLessTc,cLessTc]=dot(0.5,(Ap - dot(Ap,cos(dot(pi / (tc),(t(n) - t0s(rLessTc,cLessTc)))))))
# ../decompositionproject_python/activation_lib\activationDynamics.m:194
                            neuralMatrix[rInTcThr,cInTcThr]=dot(0.5,(Ap - dot(Ap,cos(dot(pi / (dot(2,(thr - tc))),(t(n) - t0s(rInTcThr,cInTcThr) + dot(2,(thr - dot(3,tc) / 2))))))))
# ../decompositionproject_python/activation_lib\activationDynamics.m:195
                            neuralActivation[n,arange()]=sum(neuralMatrix(arange(1,n_spike),arange()),1,'omitnan')
# ../decompositionproject_python/activation_lib\activationDynamics.m:198
            if tetanusModel:
                neuralActivationTransf[n,arange()]=dot(normSat,(1 - dot(1,exp(dot(- tetParam,(neuralActivation(n,arange()))))))) / (1 + dot(1,exp(dot(- tetParam,(neuralActivation(n,arange()))))))
# ../decompositionproject_python/activation_lib\activationDynamics.m:201
    
    #     size(excitation,2)
#                     totalMUs*recCurve(max(force))
#                     totalMUs*recCurve(max(force))/size(excitation,2)
    # hold on, plot(t,neuralActivation)
    if tetanusModel:
        neuralActivation=copy(neuralActivationTransf)
# ../decompositionproject_python/activation_lib\activationDynamics.m:210
    
    return neuralActivation,normSat
    
if __name__ == '__main__':
    pass
    