# Generated with SMOP  0.41-beta
from libsmop import *
# ../decompositionproject_python/activation_lib\getTwitchProperties.m

    
@function
def getTwitchProperties(dischargeRates=None,recruitmentThreshold=None,muscle=None,transtable=None,*args,**kwargs):
    varargin = getTwitchProperties.varargin
    nargin = getTwitchProperties.nargin

    # getTwitchProperties estimates twitch properties based on previously found
# mapping (Gogeascoechea et al. 2023)
# tcs are contraction times and aps are peak amplitudes
# Inputs:
# - dischargeRates: Vector of mean discharge rates
# - recruitmentThreshold: Vector of mean recruitment Thresholds
# - muscle: string -> 'TA' 'SOL' 'GASlat' 'GASmed' 'PERLONGUS' 'PERTERT'
# - transtable: load from translationTable.mat
    
    # preparing data to be in columns
    if size(dischargeRates,2) > size(dischargeRates,1):
        dischargeRates=dischargeRates.T
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:13
    
    if size(recruitmentThreshold,2) > size(recruitmentThreshold,1):
        recruitmentThreshold=recruitmentThreshold.T
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:16
    
    # find id for desired muscle
    m=find(strcmp(muscle,transtable.muscles))
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:20
    # find coefficients on table
    for n in arange(1,length(transtable.columnsTransTable)).reshape(-1):
        param=transtable.columnsTransTable[n]
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:24
        if 'dr0' == param:
            dr0=transtable.translationTable(m,n)
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:27
        else:
            if 'rt0' == param:
                rt0=transtable.translationTable(m,n)
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:29
            else:
                if 'c1' == param:
                    c1=transtable.translationTable(m,n)
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:31
                else:
                    if 'c2' == param:
                        c2=transtable.translationTable(m,n)
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:33
                    else:
                        if 'p_tc' == param:
                            p_tc=transtable.translationTable(m,arange(n,n + 1))
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:35
                        else:
                            if 'p_ap' == param:
                                p_ap=transtable.translationTable(m,arange(n + 1,n + 2))
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:37
    
    coeff=concat([[c1],[c2]])
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:41
    tcs=dot(dot(p_tc(1),concat([dischargeRates / 40 - dr0,recruitmentThreshold - rt0])),coeff) + p_tc(2)
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:43
    aps=dot(dot(p_ap(1),concat([dischargeRates / 40 - dr0,recruitmentThreshold - rt0])),coeff) + p_ap(2)
# ../decompositionproject_python/activation_lib\getTwitchProperties.m:44
    return tcs,aps
    
if __name__ == '__main__':
    pass
    