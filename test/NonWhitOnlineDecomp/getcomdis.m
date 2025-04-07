function [comdis, roa, coordinates_com] = getcomdis(PulseT,PulseT2, distime, distime2, maxlag, jitter, tol, fsamp)

jit = round(jitter*fsamp);
% Generate binary spike trains
firings = zeros(size(PulseT));

for i = 1:size(PulseT,1)
    firings(i , distime{i}) = 1;
    distimmp{i} = [];
    for j = 1:jit
        distimmp{i} = [distimmp{i} distime{i}-j];
        distimmp{i} = [distimmp{i} distime{i}+j];
    end
    distimmp{i} = [distimmp{i} distime{i}];
end
firings = firings(:, 1:length(PulseT));

firings2 = zeros(size(PulseT2));
for i = 1:size(PulseT2,1)
    firings2(i , distime2{i}) = 1;
    distimmp2{i} = [];
    for j = 1:jit
        distimmp2{i} = [distimmp2{i} distime2{i}-j];
        distimmp2{i} = [distimmp2{i} distime2{i}+j];
    end
    distimmp2{i} = [distimmp2{i} distime2{i}];
end
firings2 = firings2(:, 1:length(PulseT2));

comdis = zeros(size(PulseT,1), size(PulseT2,1));
roa = zeros(size(PulseT,1), size(PulseT2,1));
num_com = 0;
% Loop through each combination of spike trains
for i = 1:size(PulseT, 1)
    for j = 1:size(PulseT2, 1)

        % Cross-correlation between spike trains
        [c, lags] = xcorr(firings(i, :), firings2(j, :), maxlag * 2, 'normalized');
        [correl, idx] = max(c);
        
        % Determine the common discharge
        if correl > 0.2
            distimetemp = distimmp2{j} + lags(idx);
        else
            distimetemp = distimmp2{j};
        end
        
        com = intersect(distimmp{i}, distimetemp);
        com([false, diff(com) == 1]) = [];
        comdis(i, j) = length(com) / max([length(distime{i}), length(distime2{j})]);

        if comdis(i, j) > 0.3
            coordinates_com(num_com+1, :) = [i, j]; % Track coordinates of common discharges
            Aj = length(com);
            Ij = length(distime{i}) - Aj;
            Sj = length(distime2{j}) - Aj;
            roa(i, j) = Aj / (Aj + Ij + Sj);
            RoA_com(num_com+1) = Aj / (Aj + Ij + Sj);
            num_com = num_com + 1;
        end
    end
end

end
