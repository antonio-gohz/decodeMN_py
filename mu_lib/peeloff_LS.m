function [EMGres,What] = peeloff_LS(spikeTrains, EMG, nw)
%  [What,wwsigs]=estimWaveforms(X, Y, nw)
%
%  Computes estimate of spike waveform for cells from spike trains and
%  electrode data, using (correct) least-squares regression
%
%  Input:  
%  ------
%   spikeTrains [nsamps x ncells] - each column holds spike train of single neuron
%   EMG [nsamps x nelec] -  electrode data
%   nw [1 x 1] - number of time bins in the spike waveform
%
%  Output:  
%  -------
%   What [nw x ne x ncells] - estimated waveforms for each cell
%   sigs [ncells x 1] - posterior stdev of each neuron's waveform coefficients
%
% Adapated from: jw pillow 8/18/2014 --> BinaryPursuitSpikeSorting-master
transposeEMG =0;
if size(EMG,2)>size(EMG,1)
    EMG=EMG';
    transposeEMG =1;
end

[nt,nc] = size(spikeTrains); % number of time bins and number of cells
ne = size(EMG,2); % number of electrodes
nw2 = nw/2;

% Compute blocks for covariance matrix XX and cross-covariance XY
XXblocks = zeros(nc*nw,nc);
XY = zeros(nc*nw,ne);
X_toeplitz =  zeros(nt,nc*nw);
X_toeplitz(:,1:nc) =  spikeTrains;

for jj = 1:nw
    inds = ((jj-1)*nc+1):(jj*nc);
    XXblocks(inds,:) = spikeTrains(1:end-jj+1,:)'*spikeTrains(jj:end,:); % spike train covariance
    XY(inds,:) = spikeTrains(max(1,nw/2-jj+2):min(nt,nt-jj+nw2+1),:)'*...
        EMG(max(1,jj-nw2):min(nt,nt+jj-nw2-1),:); % cross-covariance
    X_toeplitz(jj:end,inds) = spikeTrains(1:end-jj+1,:);%X(1:length(X_toeplitz(jj:end-jj+1,inds)),:);
end

% Insert blocks into covariance matrix XX
XX = zeros(nc*nw,nc*nw);
for jj = 1:nw
    inds1 = ((jj-1)*nc+1):(nc*nw);
    inds2 = ((jj-1)*nc+1):(jj*nc);
    XX(inds1,inds2) = XXblocks(1:(nw-jj+1)*nc,:);  % below diagonal blocks
    XX(inds2,inds1) = XXblocks(1:(nw-jj+1)*nc,:)'; % above diagonal blocks
end
What = XX\XY; % do regression

EMGrec = X_toeplitz*What;
EMGrec = [EMGrec(nw2+1:end,:);zeros(nw2,size(EMGrec,2))];
EMGres = EMG-EMGrec;
% Note: the "correct" formula should be diag(inv(Xcov)), but this is
% close, ignoring edge effects, and much faster;

if nargout > 1
    What = permute(reshape(What,[nc,nw,ne]),[3,2,1]);  % reshape into tensor [2,3,1]
end
if transposeEMG
    EMGres=EMGres';
end