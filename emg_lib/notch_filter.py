from scipy.signal import filtfilt, iirnotch

def notch_filter(data, fs, notch_freq=50, order=4, axis=-1):
    nyquist = 0.5 * fs
    freqs = [notch_freq / nyquist]
    
    b, a = iirnotch(freqs[0], order)
    return filtfilt(b, a, data, axis=axis)  # Apply along the given axis
