import numpy as np
from decomposition_lib.getspikes import getspikes
# from plots.updateLivePlots import updateLivePlots

def maximizeSIL(X, w, spikes,sil, ipt, norm_ipt, centroids, fsamp, showPlots=False, live_plot=None):
    """
    Optimization loop to maximize the silhouette of motor unit discharges.

    Parameters:
    - spikes (np.array): Initial discharge times
    - X (np.array): Whitened signal
    - sil (float): Initial silhouette value
    - w (np.array): Initial weights (MU filter)
    - fsamp (float): Sampling frequency
    - showPlots (bool): Whether to show live plots (optional)
    - h, ax: Plot handles for updating live plots (optional)

    Returns:
    - w_last (np.array): Optimized weights (MU filter)
    - spikes_last (np.array): Optimized discharge times of the motor unit
    - sil_last (float): Final silhouette value
    - ipt_last (np.array): Last computed MU pulse train
    - norm_ipt (float): Normalization factor for MU pulse train
    - centroids (np.array): Cluster centroids of the pulse train
    """ 
    k = 1
    sil_last = np.finfo(float).eps  # Smallest possible value for stability in comparison
    #ipt, spikes, sil, norm_ipt, centroids = getspikes(w, X, fsamp)

    # For plotting
    x_x = np.linspace(0, X.shape[1] / fsamp, X.shape[1])
    x_w = np.arange(len(w))

    while sil > sil_last:
        # Save last values
        sil_last = sil
        spikes_last = spikes
        w_last = w
        ipt_last = ipt

        # Optional live plotting
        # if showPlots:
        #     h, ax = updateLivePlots(
        #         [[x_x, ipt_last / np.max(ipt_last)], 
        #          [x_x[spikes], ipt_last[spikes] / np.max(ipt_last)], 
        #          [x_w, w]], h, ax)
        #     ax[0].set_title(f"2 Refinement SIL: {sil_last:.4f} Iteration: {k}")

        if showPlots and live_plot is not None:
            # Check if the plot is still open before updating
            if live_plot.is_open:
                data = [[(x_x, ipt_last / np.max(ipt_last))], \
                        [(x_x[spikes],ipt_last[spikes] / np.max(ipt_last)), (x_w,w)]]  # List of lists for the axes
                # Update each subplot's line with new data
                live_plot.update(data) 
            else:
                print("Plot closed, stopping updates.")
                live_plot.close()  # Gracefully close the plot
                live_plot = None  # Set live_plot to None if the plot is closed

        # Update weights based on the last detected spikes
        k += 1
        w = np.sum(X[:, spikes_last], axis=1)
        w = w / np.linalg.norm(w)

        # Recalculate discharge times and silhouette
        ipt, spikes, sil, norm_ipt, centroids = getspikes(w, X, fsamp)

    # Ensure discharge times are saved if the last spike count is low
    if len(spikes_last) < 2:
        _, spikes_last = getspikes(w, X, fsamp)

    return w_last, spikes_last, sil_last, ipt_last, norm_ipt, centroids, live_plot

    