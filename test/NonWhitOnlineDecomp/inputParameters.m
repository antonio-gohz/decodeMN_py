function [eSIG, parameters] = inputParameters(EMG, decompParameters, signal, whiten, tanh_denoise,fs,cut)
%   InputParameters:  Prepare data from offline decomposition to get online
%   parameters
%   [eSig, parameters] = inputParameters(EMGFILT, 'decompParameters', signal, 'whiten', tanh_true, fs)
%   Prepares the output of the offline decomposition for extraction of
%   either whitened or non-whitened online parameters.
%
%   Input:
%   - EMG: EMG observations
%   - 'decompParameters': contains MUfilters, whitening matrix and EMGmask
%   - 'signal': output signal of offline decomposition with the discharge
%   times
%   - 'whiten': flag to use whitened method in extraction of online
%   parameters
%   - 'tanh_true': flag to use tanh_denoising on extended non-whitened data
%   - 'fs': sample frequency of EMG data
%   - 'cut': flag to cut the signals to avoid transients in pulse trains
%     (After extend() the first and last #extensionfactor# columns contain
%     a lot of zeros which causes transients in the estimated pulse trains.
%     The first and last second will be cut off to avoid these transients
%     if cut > 0)
%
%   Output:
%   - eSIG: Extended EMG signals, depending on the inputs it will be
%   (non)-whitened 
%   - parameters: Struct containing decomposition parameters for online
%   biofeedback
    
% Whitened method of CKC, so covariance matrix C_zz = I
    if whiten == 1
        %Extend observations, exlude masked channels
        eSIG = extend(EMG(decompParameters.EMGmask,:),decompParameters.extensionfactor);
        %whiten the observations with the whitening matrix from offline
        %decomposition
        eSIG = decompParameters.whitMat * eSIG;
        
        
       if tanh_denoise > 0
            normtanh =  tanh_denoise*std(eSIG,[],3);
            eSIG = normtanh.*tanh(eSIG./normtanh);
        end
        %Cut of the edges of the signal, same as in offline decomposition.
        %This step is necessary to avoid a transient in the beginning and
        %end of the signal.
        if cut > 0
            eSIG = eSIG(:,round(fs):end-round(fs));
        end

        %Use whitened MUfilters
        parameters.MUfilters = decompParameters.MUFilters;
        %Use inverse covariance matrix of whitened signals C_{zz}^{-1} = I
        parameters.ReSIG = eye(size(decompParameters.dewhitMat));
        parameters.iReSIGt = eye(size(decompParameters.dewhitMat));
        % Save extensionfactor
        parameters.extensionfactor = decompParameters.extensionfactor;
    
% Non whitened method of CKC, so covariance matrix C_xx
    else 
        %Extend observations, exlude masked channels
        eSIG = extend(EMG(decompParameters.EMGmask,:),decompParameters.extensionfactor);
        
        if tanh_denoise > 0
            normtanh =  tanh_denoise*std(eSIG,[],2);
            eSIG = normtanh.*tanh(eSIG./normtanh);
        end

        %Cut of the edges of the signal, same as in offline decomposition.
        %This step is necessary to avoid a transient in the beginning and
        %end of the signal.
        if cut > 0
            eSIG = eSIG(:,round(fs):end-round(fs));
        end
        %Use non-whitened MUfilters by averaging values at discharge times
        %of extended signals eSIG
        parameters.MUFilters = getMUfilters(eSIG, signal.Distime);
        %Use the inverse covariance matrix of non-whitneed signals C_{xx}^{-1}
        parameters.ReSIG = eSIG*eSIG'/length(eSIG);
        parameters.iReSIGt = pinv(parameters.ReSIG);

       % Save extensionfactor
        parameters.extensionfactor = decompParameters.extensionfactor;
    end
end

