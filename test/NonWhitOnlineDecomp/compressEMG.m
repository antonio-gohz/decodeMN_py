function EMG_compressed = compressEMG(EMG)
    % Function to compress EMG data based on standard deviation
    
    % Initialize the compressed matrix
    EMG_compressed = EMG;
    
    % Calculate the standard deviation for each row
    std_dev = std(EMG, 0, 2); % Compute standard deviation along rows (dimension 2)
    
    % Define a scaling factor based on standard deviation
    scaling_factor = 1 ./ (1 + std_dev); % Adjust as needed
    
    % Apply hyperbolic tangent transformation to each row with dynamic scaling factor
    for i = 1:size(EMG, 1)
        % Extract the current row
        data = EMG(i, :);
        
        % Calculate attenuation factor based on standard deviation
        attenuation_factor = scaling_factor(i); % Use the scaling factor for this row
        
        % Apply hyperbolic tangent transformation
        data_tanh = tanh(data * attenuation_factor) .* data;
        
        % Store the compressed row back into the matrix
        EMG_compressed(i, :) = data_tanh;
    end
end
