% Example covariance matrix
exFactor = 18;
eSIG = extend(EMG(emgMask,:),exFactor);
eSIG = eSIG(:,1000:end-1000);
[E, D] = pcaesig(eSIG);
[wSIG, whiteningMatrix, dewhiteningMatrix] = whiteesig(eSIG, E, D);

Cxx = eSIG*eSIG'/length(eSIG);

% Compute K(Czz)
Czz_calc = whiteningMatrix * Cxx * whiteningMatrix';
K_Czz = computeKLD(eye(4));

function K_Czz = computeKLD(Czz)
    % Function to compute the Kullback-Leibler divergence K(Czz)
    % Input:
    %   Czz - Covariance matrix (n x n)
    % Output:
    %   K_Czz - The computed value of K(Czz)
    
    % Check if the input is a square matrix
    [n, m] = size(Czz);
    if n ~= m
        error('Czz must be a square matrix.');
    end

    
    % Compute the trace of Czz
    trace_Czz = trace(Czz);
    
    % Compute the log determinant of Czz
    log_det_Czz = log(det(Czz));
    % Compute K(Czz)
    K_Czz = 0.5 * (trace_Czz - log_det_Czz - n);
end

