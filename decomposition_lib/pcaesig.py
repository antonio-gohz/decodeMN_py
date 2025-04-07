import numpy as np
from sklearn.decomposition import PCA

def pcaesig(signal, n_components=None):
    """
    Perform PCA on the input signal using sklearn.

    Parameters:
    signal (np.ndarray): The input signal, expected to be row-wise.
    n_components (int, optional): Number of components to keep. If None, all components are kept.

    Returns:
    E (np.ndarray): The matrix whose columns are the corresponding eigenvectors.
    D (np.ndarray): The diagonal matrix of eigenvalues.
    """
    pca = PCA(n_components=n_components, whiten=False)
    pca.fit(signal.T)  # Transpose to fit PCA

    # Eigenvalues and eigenvectors
    E = pca.components_.T  # Shape (n_features, n_components)
    D = np.diag(pca.explained_variance_)  # Shape (n_components, n_components)

    # Compute rank tolerance
    eigenvalues = np.diag(D)
    rankTolerance = np.mean(eigenvalues[len(eigenvalues) // 2:]) if len(eigenvalues) > 1 else 0

    # Set rank tolerance to zero if negative
    if rankTolerance < 0:
        rankTolerance = 0

    # Select columns based on the rank tolerance
    mask = eigenvalues > rankTolerance
    E = E[:, mask]  # Keep only relevant eigenvectors
    D = D[mask][:, mask]  # Keep only relevant eigenvalues in diagonal form


    return E, D