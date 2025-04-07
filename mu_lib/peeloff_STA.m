%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peel off of the Motor unit spike train

% Input: 
%   X = whitened signal
%   spikes = discharge times of the motor unit
%   fsamp = sampling frequency
%   win = window to identify the motor unit action potential with spike trigger averaging


% Output:
%   X = residual of the whitened signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = peeloff_STA(X, firings, fsamp, win)
transposeEMG=0;
if size(X,2)<size(X,1)
    X=X';
    transposeEMG =1;
end

windowl = round(win*fsamp);
waveform = zeros(1,windowl*2+1);

EMGtemp = zeros(size(X));
for l = 1:size(X,1)
    temp = cutMUAP(find(firings),windowl,X(l,:));
    waveform = mean(temp,1);
    EMGtemp(l,:) = conv(firings,waveform,'same');
end 
    
X = X - EMGtemp;
if transposeEMG
    X = X';
end
