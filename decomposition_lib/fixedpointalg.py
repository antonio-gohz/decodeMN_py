import numpy as np
# import matplotlib.pyplot as plt
# from plots.live_plots import LivePlot


def fixedpointalg(w, X, B, maxiter, contrastfunc, fsamp, showPlots=False, live_plot=None):
    k = 0
    delta = np.ones(maxiter)
    TOL = 0.0001  # tolerance between two iterations
    BBT = B @ B.T
    
    # Define the contrast function based on the given type
    if contrastfunc == 'square':
        gp = lambda x: x
        g = lambda x: (x**2) / 2
    elif contrastfunc == 'skew':
        gp = lambda x: (2 * x**2) / 3
        g = lambda x: (x**3) / 3
    elif contrastfunc == 'logcosh':
        gp = lambda x: np.tanh(x)
        g = lambda x: np.log(np.cosh(x))
    elif contrastfunc == 'logcosh_mod':
        gp = lambda x: np.tanh(2 * x)
        g = lambda x: np.log(np.cosh(2 * x)) / 2
    elif contrastfunc == 'pow3':
        gp = lambda x: x**3
        g = lambda x: (x**4 - 3) / 4

    # For plotting
    x_x = np.linspace(0, X.shape[1] / fsamp, X.shape[1])
    x_w = np.arange(len(w))


    while delta[k] > TOL and k < maxiter:
        # Update weights
        wlast = w.copy()  # Save last weights

        # Contrast function
        wTX = w.T @ X
        A = np.mean(gp(wTX))
        w = X @ (g(wTX)) / len(wTX) - A * w  # Update weights

        # Orthogonalization
        w = w - BBT @ w

        # Normalization
        w /= np.linalg.norm(w)

        # Update convergence criteria
        k += 1
        if k < maxiter:  # Check if k is still within bounds
            delta[k] = abs(w.T @ wlast - 1)

        # If plots are enabled, update them with current data
        if showPlots and live_plot is not None:
            # Check if the plot is still open before updating
            if live_plot.is_open:
                data = [[(x_x, wTX / np.max(wTX))], \
                        [(x_w, np.nan * np.ones_like(w)), (x_w, w)]]  # List of lists for the axes
                # Update each subplot's line with new data
                live_plot.update(data)
            else:
                print("Plot closed, stopping updates.")
                live_plot.close()  # Gracefully close the plot
                live_plot = None  # Set live_plot to None if the plot is closed

    return w, live_plot

    