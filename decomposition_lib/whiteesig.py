import numpy as np

def whiteesig(signal, E, D,  pad_length=0):
    """
    Apply whitening to the input signal.

    Parameters:
    signal (np.ndarray): The input signal to be whitened.
    E (np.ndarray): The matrix of eigenvectors from PCA.
    D (np.ndarray): The diagonal matrix of eigenvalues from PCA.

    Returns:
    whitensignals (np.ndarray): The whitened signals.
    whiteningMatrix (np.ndarray): The whitening transformation matrix.
    dewhiteningMatrix (np.ndarray): The de-whitening transformation matrix.
    """
    # Apply zero-padding if pad_length is greater than 0
    if pad_length > 0:
        original_shape = signal.shape
        signal = np.pad(signal, ((0, 0), (pad_length, pad_length)), mode='constant')

    # Calculate the whitening transformation matrix
    whiteningMatrix = E @ np.linalg.inv(np.sqrt(D)) @ E.T
    dewhiteningMatrix = E @ np.sqrt(D) @ E.T
    whitensignals = whiteningMatrix @ signal

    # Remove padding to return the whitened signal to its original length
    if pad_length > 0:
        whitensignals = whitensignals[:, pad_length:-pad_length]

    return whitensignals, whiteningMatrix, dewhiteningMatrix