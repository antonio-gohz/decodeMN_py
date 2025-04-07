import numpy as np
from scipy.signal import butter, filtfilt, iirnotch
from emg_lib import butter_filter
from emg_lib.notch_filter import notch_filter

def filter_emg(emg, fs, order=4, bandpass_lowcut=10, bandpass_highcut=500, 
               notch_freq=50, notch_harmonics=[50, 100, 150, 200, 250, 300, 350, 400, 450, 500], 
               check_for_noise=True, dimension=None):
    """
    Function to filter EMG signals, check for noisy channels, and return the EMG envelope.
    
    Parameters:
    - emg: numpy array of shape (n_channels, n_samples) representing the EMG signal.
    - fs: Sampling frequency in Hz.
    - order: Order of the filter (default: 4).
    - bandpass_lowcut: Lower cut-off frequency for bandpass filter (default: 10 Hz).
    - bandpass_highcut: Upper cut-off frequency for bandpass filter (default: 500 Hz).
    - notch_freq: Center frequency for the notch filter (default: 50 Hz).
    - notch_harmonics: List of harmonics for the notch filter (default: [50, 100, 150, ..., 500 Hz]).
    - check_for_noise: Flag to check if the input EMG has noisy channels (default: True).
    - dimension: Axis along which to apply the filter ('auto' will decide based on input matrix shape).
    
    Returns:
    - filtered_emg: The filtered EMG signal.
    - emg_mask: A boolean mask indicating good channels (True for good, False for noisy).
    - emg_envelope: The EMG envelope, which is the low-pass filtered version of rectified filtered EMG.
    """
    
    # Determine the axis based on the dimension parameter or automatically based on shape
    if dimension is None:
        # Automatically detect which dimension is larger (time vs. channels)
        dimension = 0 if emg.shape[0] > emg.shape[1] else 1
    
    # Apply bandpass filter (default 10-500 Hz) to all channels (along the determined dimension)
    filtered_emg = butter_filter(emg, bandpass_lowcut, bandpass_highcut, fs, order, 'bandpass', axis=dimension)
    
    # Apply notch filter for the given harmonics (default 50 Hz and harmonics)
    for harmonic in notch_harmonics:
        filtered_emg = notch_filter(filtered_emg, fs, harmonic, order, axis=dimension)
    
    # Check for noisy channels: Example approach (you can customize this)
    # For example, if the signal has a high standard deviation (noisy signal)
    emg_mask = np.std(filtered_emg, axis=1) > 0.1  # Threshold for detecting noisy channels
    
    # Calculate the EMG envelope (low-pass filter of the rectified signal)
    rectified_emg = np.abs(filtered_emg)  # Rectification
    emg_envelope = butter_filter(rectified_emg, 0, 10, fs, order, 'lowpass', axis=dimension)  # Low-pass filter
    
    return filtered_emg, emg_mask, emg_envelope

# Example usage:
# emg: input EMG data (n_channels, n_samples)
# fs: sampling frequency (e.g., 1000 Hz)
# filtered_emg, emg_mask, emg_envelope = filter_emg(emg, fs)