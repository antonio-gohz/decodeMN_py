%% tanh before whitening 
% Signal Extension
exFactor = 18;
eSIG = extend(EMGtrain(emgMask,1:20481),exFactor);
tanh_denoise = [1.5];

if ~isempty(tanh_denoise)
    normtanh =  tanh_denoise*std(eSIG,[],2);
    % if optimizeTanh
    % [normtanh, eSIG] = denoiseTanH(eSIG,1);
    % end
eSIGbefore= normtanh.*tanh(eSIG./normtanh);
end

% Signal Whitening
[E, D] = pcaesig(eSIG);
[wSIGbefore, whiteningMatrix, dewhiteningMatrix] = whiteesig(eSIGbefore, E, D);

eSIG = extend(EMGtrain(emgMask,1:20481),exFactor);
tanh_denoise = [1.5];

% Signal Whitening
[E, D] = pcaesig(eSIG);
[wSIGafter1, whiteningMatrix, dewhiteningMatrix] = whiteesig(eSIG, E, D);

if ~isempty(tanh_denoise)
    normtanh =  tanh_denoise*std(wSIGafter1,[],2);
    % if optimizeTanh
    % [normtanh, wSIGafter] = denoiseTanH(wSIGafter,1);
    % end
wSIGafter= normtanh.*tanh(wSIGafter1./normtanh);
end
sgtitle('Preprocessing of EMG(1,:)');
%% Covariance matrix 
CovMatrix = wSIGafter1*wSIGafter1'/length(wSIGafter1);
%% Other types of denoising before whitening
for i = 1:size(EMG, 1)
    % Extract the current row
    data = EMG(i, :);
    
    % Calculate Q1 and Q3
    Q1 = prctile(data, 5);
    Q3 = prctile(data, 95);
    
    % Calculate IQR
    IQR = Q3 - Q1;
    
    % Define lower and upper bounds
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
    
    % Replace outliers with the nearest bound
    data_iqr = data;
    data_iqr(data < lower_bound) = lower_bound;
    data_iqr(data > upper_bound) = upper_bound;
    
    % Store the attenuated row back into the matrix
    EMG_iqr(i, :) = data_iqr;
end

% Plot the first row of original and attenuated data
figure;
plot(EMG(1, :), 'DisplayName', 'Original Data');
hold on;
plot(EMG_iqr(1, :),'DisplayName', 'Attenuated Data');
hold off;

% Add titles and labels
title('EMG Data - First Row');
xlabel('Sample Index');
ylabel('Amplitude');
legend show;
grid on;
%% Compare
fs= 2048;
time = 0:1/fs:(length(EMGtrain(emgMask,1:20481))-1)./fs ;
subplot(4, 1, 1);
plot(time, EMGtrain(1, 1:20481));
title('Original data');
xlabel('Time');
ylabel('Amplitude');

% Plot wSIGbefore
subplot(4, 1, 2);
plot(time, wSIGbefore(1, 1:end-17));
title('Tanh Denoising before Whitening');
xlabel('Time');
ylabel('Amplitude');

% Plot wSIGafter
subplot(4, 1, 3);
plot(time, wSIGafter1(1, 1:end-17));
title('No denoising only whitening');
xlabel('Time');
ylabel('Amplitude');

% Plot wSIGafter
subplot(4, 1, 4);
plot(time, wSIGafter(1, 1:end-17));
title('Tanh Denoising after Whitening');
xlabel('Time');
ylabel('Amplitude');
% Display the plots
disp('Displaying the first channels of wSIGbefore and wSIGafter');
sgtitle('Preprocessing of EMG(1,:)');
%% 
correl= xcorr(eSIGbefore(1, 1:end-17),eSIG(1, 1:end-17), 'normalized');
plot(correl)
max(correl)
%% correlation analysis
tanh_denoise = 3;
normtanh =  tanh_denoise*std(eSIG,[],2);
eSIGafter= normtanh.*tanh(eSIG./normtanh);

[R,PValue] = corrplot([eSIG(1,:); eSIG(1,:)]');


