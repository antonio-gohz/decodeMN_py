import numpy as np
from scipy.signal import butter, filtfilt

def butter_filter(data, lowcut, highcut, fs, order=4, filter_type='bandpass', axis=-1):
    nyquist = 0.5 * fs
    low = lowcut / nyquist
    high = highcut / nyquist

    # Bandpass filter
    if filter_type == 'bandpass':
        b, a = butter(order, [low, high], btype='band')
    # Lowpass filter
    elif filter_type == 'lowpass':
        b, a = butter(order, high, btype='low')
    # Highpass filter
    elif filter_type == 'highpass':
        b, a = butter(order, low, btype='high')

    return filtfilt(b, a, data, axis=axis)  # Apply along the given axis
