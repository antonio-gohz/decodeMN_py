import numpy as np
from scipy.signal import correlate
from tqdm import tqdm  # For progress bar

def remduplicates(PulseT, distime, distime2, maxlag, jitter, tol, fsamp, sil=None):
    """
    Removes duplicate motor units based on shared discharge times.

    Parameters:
        PulseT (np.ndarray): Pulse train of each MU.
        distime (list): Discharge times of the motor units.
        distime2 (list): Realigned discharge times of the motor units.
        maxlag (int): Maximum lag between motor unit spike trains.
        jitter (float): Tolerance in seconds for the estimation of discharge times.
        tol (float): Threshold of shared discharge times to define duplicates.
        fsamp (float): Sampling frequency.
        sil (list, optional): Silhouette values.

    Returns:
        Pulsenew (np.ndarray): Pulse train of non-duplicated MU.
        distimenew (list): Discharge times of non-duplicated MU.
        idsnewVec (list): IDs of non-duplicated MUs.
    """
    jit = int(round(jitter * fsamp))
    firings = np.zeros_like(PulseT)

    # Create jittered discharge times
    distimmp = []
    for i in range(PulseT.shape[0]):
        firings[i, distime2[i]] = 1
        jittered_times = set(distime2[i])
        for j in range(1, jit + 1):
            jittered_times.update(distime2[i] - j)
            jittered_times.update(distime2[i] + j)
        distimmp.append(sorted(jittered_times))
    
    MUn = len(distime2)
    idsnewVec, distimenew, Pulsenew = [], [], []
    idsnew = list(range(MUn))

    for i in tqdm(range(MUn), desc="Removing duplicates"):
        if not distimmp:
            break
        
        comdis = np.zeros(len(distimmp))
        for j in range(1, len(distimmp)-1):
            corr = correlate(firings[0, :], firings[j, :], mode='valid')
            correl_max = np.max(corr)
            lag_idx = np.argmax(corr) - maxlag
            distimetemp = np.array(distimmp[j]) + lag_idx if correl_max > 0.2 else distimmp[j]
            common = np.intersect1d(distimmp[0], distimetemp)
            comdis[j] = len(common) / max(len(distime[0]), len(distime[j]))

        duplicates = [0] + [j for j in range(1, len(comdis)) if comdis[j] >= tol]

        if len(duplicates) > 1:
            active_period = np.array([(distime[d][-1] - distime[d][0]) / fsamp for d in duplicates])
            ranks1 = np.argsort(np.argsort(active_period))
            
            if sil is None:
                isi_metric = [np.std(np.diff(distime[d])) / np.mean(np.diff(distime[d])) for d in duplicates]
                ranks2 = np.argsort(np.argsort(isi_metric))
            else:
                ranks2 = np.argsort(sil[duplicates])
                sil = [sil[d] for d in range(len(sil)) if d not in duplicates]

            ranks = ranks1 * np.mean(np.diff(active_period)) / np.mean(active_period) + \
                    ranks2 * np.mean(np.diff(isi_metric)) / np.mean(isi_metric)
            survivor = duplicates[np.argmax(ranks)]
        else:
            survivor = duplicates[0]

        idsnewVec.append(idsnew[survivor])
        distimenew.append(distime[survivor])
        Pulsenew.append(PulseT[survivor, :])

        # Remove duplicates
        for idx in sorted(duplicates, reverse=True):
            del distime[idx]
            del distime2[idx]
            del distimmp[idx]
            del idsnew[idx]
            PulseT = np.delete(PulseT, idx, axis=0)
            firings = np.delete(firings, idx, axis=0)

    return np.array(Pulsenew), distimenew, idsnewVec
