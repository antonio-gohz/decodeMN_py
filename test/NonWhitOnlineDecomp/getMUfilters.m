function MUFilters = getMUfilters(eSIG, Distime)
%   MUfilters: Extracts MUfilters based on extended observations (eSIG) and
%   dischargetimes of motor units 
% 
%   Input:
%   'eSIG': extended non-whitened observations 
%   'Distime': 1 x nMU cell which contains discharge times of different motor units
%
%   Output:
%   'MUfilters': mean value of eSIG at discharge times (%based on online
%   method of Yeung et al.)  

% Recalculate MUfilters
MUFilters = zeros(size(eSIG,1), length(Distime));
for i = 1:length(Distime)
    MUFilters(:,i) =  mean(eSIG(:,Distime{i}),2);
end