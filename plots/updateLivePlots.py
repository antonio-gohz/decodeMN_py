import matplotlib.pyplot as plt
import numpy as np

def updateLivePlots(data, lines, ax, xlimits=None):
    for i, line_set in enumerate(lines):
        for j, line in enumerate(line_set):
            if data[i][1] is not None:
                line.set_xdata(data[i][0])
                line.set_ydata(data[i][1][j])
    
    # Update axis limits if provided
    if xlimits is not None:
        for a in ax:
            a.set_xlim(xlimits)
    
    plt.draw()
    plt.pause(0.001)  # Small pause to allow the figure to update
