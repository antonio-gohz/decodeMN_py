import numpy as np
from decomposition_lib.getspikes import getspikes
from plots.updateLivePlots import updateLivePlots

def minimizeCOVISI(w, X, cov, fsamp, show_plots=False, h=None, ax=None):
    """
    Optimization loop to minimize the CoV of interspike intervals for MU filter weights.
    
    Parameters:
    - w (np.array): Initial weights
    - X (np.array): Whitened signal
    - cov (float): Initial CoV of interspike intervals
    - fsamp (float): Sampling frequency
    - show_plots (bool): Whether to show live plots (optional)
    - h, ax: Plot handles for updating live plots (optional)
    
    Returns:
    - w_last (np.array): Updated weights (MU filter)
    - spikes_last (np.array): Discharge times of the motor unit
    - cov_last (float): Updated CoV of the interspike intervals
    - sil_last (float): Silhouette value from the last iteration
    - ipt_last (np.array): Last computed MU pulse train
    - norm_ipt (float): Normalization factor for MU pulse train
    - centroids (np.array): Cluster centroids of the pulse train
    """
    k = 1
    cov_last = cov + 0.1
    ipt, spikes, sil, norm_ipt, centroids = getspikes(w, X, fsamp)
    
    # For plotting
    x_x = np.linspace(0, len(X) / fsamp, len(X))
    x_w = np.arange(len(w))
    
    while cov < cov_last:
        # Save last values
        cov_last = cov
        spikes_last = spikes
        w_last = w
        sil_last = sil
        ipt_last = ipt

        # Optional live plotting
        if show_plots:
            h, ax = updateLivePlots(
                [[x_x, ipt_last / np.max(ipt_last)], 
                 [x_x[spikes], ipt_last[spikes] / np.max(ipt_last)], 
                 [x_w, w]], h, ax)
            ax[0].set_title(f"2 Refinement CoV: {cov_last:.4f} Iteration: {k}")

        # Update weights based on the last detected spikes
        k += 1
        w = np.sum(X[:, spikes_last], axis=1)
        w = w / np.linalg.norm(w)

        # Recalculate discharge times and CoV
        ipt, spikes, sil, norm_ipt, centroids = getspikes(w, X, fsamp)
        isi = np.diff(spikes / fsamp)
        cov = np.std(isi) / np.mean(isi)

    # Ensure discharge times are saved if the last spike count is low
    if len(spikes_last) < 2:
        _, spikes_last = getspikes(w, X, fsamp)
    
    return w_last, spikes_last, cov_last, sil_last, ipt_last, norm_ipt, centroids

    