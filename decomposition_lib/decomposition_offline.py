import numpy as np
from tqdm import tqdm  # For progress bar
from decomposition_lib import (
    pcaesig,
    whiteesig,
    extend,
    fixedpointalg,
    maximizeSIL,
    minimizeCOVISI,
    getspikes,
    remduplicates
)
from plots.prepareLivePlots import prepareLivePlots
from plots.updateLivePlots import updateLivePlots
from plots.live_plots import LivePlot
import matplotlib.pyplot as plt



def decomposition_offline(emg_filt, fs, exfactor,nbIterations, preOptFilters=None, refineStrategy='SIL', \
                              showPlots=False, h=None, ax=None, peeloff_flag=False, removeDuplicates=True, qc_threshold=0.85):
    """
    FastICA implementation with refinement steps to extract independent components.

    Parameters:
    X (np.ndarray): Input data matrix (features x samples).
    nbIterations (int): Number of iterations for the FastICA algorithm.
    preOptFilters (np.ndarray or None): Pre-optimized filters for initialization, if available.
    refineStrategy (str): Strategy for refining the output ('SIL' or 'COV').
    fs (float): Sampling frequency.
    showPlots (bool): Flag to show plots during processing.
    h, ax: Plotting handles, if needed.
    peeloff_flag (bool): Flag to indicate peeling off reliable sources.
    removeDuplicates (bool): Flag to indicate if duplicates should be removed.

    Returns:
    MUFilters (np.ndarray): Matrix of independent components.
    SIL (np.ndarray): Signal-to-Interference Ratio for each component.
    COV (np.ndarray): Coefficient of Variation for each component.
    PulseT (np.ndarray): Time of detected pulses.
    Distime (list): List of detected time indices.
    normIPT (np.ndarray): Normalization IPT values.
    centroids (np.ndarray): Centroid positions for each component.
    """

    esig = extend(emg_filt, exfactor)
    # Whiten the signal using PCA
    E, D = pcaesig(esig)
    wesig, whitening_matrix, dewhitening_matrix = whiteesig(esig, E, D)
    # Removing the edges
    wesig = wesig[:,exfactor:-exfactor]


    # plt.figure
    # plt.plot(wesig[0,:])
    # plt.show()
    # plt.pause(0.1)

    # Initialize flags and variables
    if preOptFilters is None:
        flagPreOptFilters = False
        Xtmp = np.copy(wesig)
        # Define percentiles for outlier thresholds
        lower_bound, upper_bound = np.percentile(Xtmp, [1, 99])

        # Replace outliers (values outside the 1st and 99th percentiles) with 0
        Xtmp[(Xtmp < lower_bound) | (Xtmp > upper_bound)] = 0        
        actind = np.sum(np.abs(Xtmp)**2, axis=0)  # Activity index
        idx1 = np.zeros(nbIterations, dtype=int)
        normIPT = np.zeros(nbIterations)
    else:
        flagPreOptFilters = True
        nbIterations = preOptFilters.shape[1]  # Number of iterations based on preOptFilters

    MUFilters = np.zeros((wesig.shape[0], nbIterations))  # Store reliable vectors
    SIL = np.zeros(nbIterations)
    COV = np.zeros(nbIterations)
    qc_metric = np.zeros(nbIterations)  # Quality control metric
    PulseT = np.zeros((nbIterations, wesig.shape[1]))
    Distime = [None] * nbIterations  # Initialize as list of None
    centroids = np.zeros((nbIterations, 2))
    # Initialize plots

    # Check if plotting is enabled and prepare plots if so
    if showPlots:
            # Initialize your live plots
            num_axes = 2  # or however many axes you need
            num_lines = [1,2]  # default is 1 line per axis
            live_plot = LivePlot(num_axes, num_lines)    
    else:
        live_plot = None

    # Progress bar for iterations
    for j in tqdm(range(nbIterations), desc='Decomposition'):
        if not flagPreOptFilters:
            idx1[j] = np.argmax(actind)
            w = wesig[:, idx1[j]]  # Initialize w
            w -= MUFilters @ MUFilters.T @ w  # Orthogonalization
            w /= np.linalg.norm(w)  # Normalization
            # Fixed point algorithm
            w, live_plot = fixedpointalg(w, wesig, MUFilters, maxiter=100, contrastfunc='logcosh', fsamp=fs, showPlots=showPlots, live_plot=live_plot)
        else:
            w = preOptFilters[:, j]  # Use pre-optimized filters for initialization

        # Maximization of SIL
        PulseT[j, :], Distime[j], SIL[j], normIPT[j], centroids[j, :] = getspikes(w, wesig, fs)

        if len(Distime[j]) > 10:  # Ensure sufficient peaks detected
            if refineStrategy == 'SIL':
                w, Distime[j], SIL[j], PulseT[j, :], normIPT[j], centroids[j, :] , live_plot = \
                    maximizeSIL(wesig, w, Distime[j], SIL[j], PulseT[j, :], normIPT[j], centroids[j, :], fs, showPlots=showPlots, live_plot=live_plot)
                ISI = np.diff(Distime[j] / fs)  # Interspike interval
                COV[j] = np.std(ISI) / np.mean(ISI)  # Coefficient of variation
                qc_metric[j] = SIL[j]
 
            elif refineStrategy == 'COV':
                ISI = np.diff(Distime[j] / fs)  # Interspike interval
                COV[j] = np.std(ISI) / np.mean(ISI)  # Coefficient of variation
                w, Distime[j], COV[j], SIL[j], PulseT[j, :], normIPT[j], centroids[j, :] , live_plot= \
                    minimizeCOVISI(w, wesig, COV[j], fs, showPlots=showPlots, live_plot=live_plot)
                qc_metric[j] = COV[j]

            # if peeloff_flag:  # Peel-off the reliable source
            #     wesig = peeloff(wesig, spikes, round(float(app.FrequencyDropDown.Value)), 0.025)

        MUFilters[:, j] = w  # Store the current filter

        actind[idx1[j]] = 0  # Remove the previous vector

    # Quality control metrics
    qc_sign = 1 if refineStrategy == 'SIL' else -1

    idsnew = np.arange(nbIterations)
    idsnew = idsnew[qc_sign * qc_metric >= qc_sign * qc_threshold]
    MUFilters = MUFilters[:, qc_sign * qc_metric >= qc_sign * qc_threshold]
    COV = COV[qc_sign * qc_metric >= qc_sign * qc_threshold]
    PulseT = PulseT[qc_sign * qc_metric >= qc_sign * qc_threshold, :]
    Distime = [d for i, d in enumerate(Distime) if qc_sign * qc_metric[i] >= qc_sign * qc_threshold]
    normIPT = normIPT[qc_sign * qc_metric >= qc_sign * qc_threshold]
    centroids = centroids[qc_sign * qc_metric >= qc_sign * qc_threshold, :]
    SIL = SIL[qc_sign * qc_metric >= qc_sign * qc_threshold]

    # Remove duplicates
    if removeDuplicates:
        PulseT, Distime, idsnewVec = remduplicates(PulseT, Distime, Distime, round(fs / 40), 0.00025, 0.3, fs)
        SIL = SIL[idsnewVec]
        COV = COV[idsnewVec]
        MUFilters = MUFilters[:, idsnewVec]
        normIPT = normIPT[idsnewVec]
        centroids = centroids[idsnewVec, :]
        idsnew = idsnew[idsnewVec]

    return MUFilters, SIL, COV, PulseT, Distime, normIPT, centroids
