function [tcs,aps] = getTwitchProperties(dischargeRates,recruitmentThreshold,muscle,transtable)
% getTwitchProperties estimates twitch properties based on previously found
% mapping (Gogeascoechea et al. 2023)
% tcs are contraction times and aps are peak amplitudes
% Inputs:
% - dischargeRates: Vector of mean discharge rates
% - recruitmentThreshold: Vector of mean recruitment Thresholds
% - muscle: string -> 'TA' 'SOL' 'GASlat' 'GASmed' 'PERLONGUS' 'PERTERT'
% - transtable: load from translationTable.mat

% preparing data to be in columns 
if size(dischargeRates,2)>size(dischargeRates,1)
    dischargeRates = dischargeRates';
end
if size(recruitmentThreshold,2)>size(recruitmentThreshold,1)
    recruitmentThreshold = recruitmentThreshold';
end

% find id for desired muscle
m = find(strcmp(muscle,transtable.muscles));

% find coefficients on table 
for n=1:length(transtable.columnsTransTable)
    param = transtable.columnsTransTable{n};
    switch param
        case 'dr0'
            dr0 = transtable.translationTable(m,n);
        case 'rt0'
            rt0 = transtable.translationTable(m,n);
        case 'c1'
            c1 = transtable.translationTable(m,n);
        case 'c2'
            c2 = transtable.translationTable(m,n);
        case 'p_tc'
            p_tc = transtable.translationTable(m,n:n+1);
        case 'p_ap'
            p_ap = transtable.translationTable(m,n+1:n+2);
    end
end

coeff = [c1;c2];
      
tcs = p_tc(1)*[dischargeRates./40-dr0, recruitmentThreshold-rt0]* coeff + p_tc(2);
aps = p_ap(1)*[dischargeRates./40-dr0, recruitmentThreshold-rt0]* coeff + p_ap(2);

end

