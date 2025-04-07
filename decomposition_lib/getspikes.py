import numpy as np
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def getspikes(w, X, fs, iReSIGt=None):
    if iReSIGt is None:  # offline case where X is whitened
        iReSIGt = np.eye(X.shape[0])

    ipt = w.T @ iReSIGt @ X
    ipt = ipt * np.abs(ipt)  # 4a: Estimate the source (abs)

    # 4b: Peak detection
    peaks, _ = find_peaks(ipt, distance=round(fs * 0.02))
    spikesTmp = peaks
    # Remove outliers
    ipt_spikes = ipt[spikesTmp]
    valid_spikes = ipt_spikes[(ipt_spikes >= np.percentile(ipt_spikes, 1)) & (ipt_spikes <= np.percentile(ipt_spikes, 99))]
    normIPT = np.mean(np.partition(valid_spikes, -10)[-10:])  # Mean of the top 10 values

    ipt = np.tanh(ipt / normIPT)

    if len(spikesTmp) > 1:
        # 4c: KMeans classification
        kmeans = KMeans(n_clusters=2)
        L = kmeans.fit_predict(ipt[spikesTmp].reshape(-1, 1))
        centroids = kmeans.cluster_centers_.flatten()
        idx2 = np.argmax(centroids)  # Spikes should be in the class with the highest centroid
        distTimes = spikesTmp[L == idx2]
        # Calculate silhouette score
        labels = kmeans.labels_
        sil = silhouette_score(ipt[spikesTmp].reshape(-1, 1), labels)

        # within = kmeans.inertia_  # Sum of squared distances to the nearest cluster center
        # # Assuming idx2 is an integer representing the index of the highest centroid
        # other_indices = np.array([0, 1])  # Indices of clusters
        # other_indices = np.setdiff1d(other_indices, idx2)  # Find the other index
        # between = np.sum(kmeans.transform(ipt[spikesTmp].reshape(-1, 1))[:, other_indices])  # Get distances to the other cluster
        # # between = np.sum(kmeans.transform(ipt[spikesTmp].reshape(-1, 1))[:, setdiff([0, 1], idx2)])  # D is not directly accessible
        # sil = (between - within) / max(within, between) if max(within, between) > 0 else 0  # Silhouette measure

    else:
        distTimes = spikesTmp
        sil = 0

    return ipt, distTimes, sil, normIPT, centroids if len(spikesTmp) > 1 else None

    