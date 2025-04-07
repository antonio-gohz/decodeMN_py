import matplotlib.pyplot as plt

def prepareLivePlots(numMUs, **kwargs):
    # Default values for optional parameters
    ylimits = kwargs.get('ylimits', [(0, 40), (0.5, numMUs + 0.5), (0, 1)])
    ylabels = kwargs.get('ylabels', ['Discharge rate (Hz)', 'MU Index', 'MU-activations (a.u.)'])
    xlabels = kwargs.get('xlabels', ['Samples'] * len(ylabels))
    xlimits = kwargs.get('xlimits', [0, 10])
    lineStyles = kwargs.get('lineStyles', ['none', '-', '-'])
    markerStyles = kwargs.get('markerStyles', ['o', 'none', 'none'])
    colors = kwargs.get('colors', [
        (0.2667, 0.4667, 0.6667),
        (0.9333, 0.4000, 0.4667),
        (0.1333, 0.5333, 0.2000),
        (0.8, 0.7333, 0.2667),
        (0.4, 0.8, 0.9333),
        (0.6667, 0.2, 0.4667),
        (0.7333, 0.7333, 0.7333)
    ])
    axesLineConfig = kwargs.get('axesLineConfig', [(0, 0), (1, 1), (2, 2)])
    linewidth = kwargs.get('linewidth', 2)

    # Create the figure and subplots
    fig, ax = plt.subplots(len(ylimits), 1, figsize=(10, 7))
    fig.suptitle('MU Feedback', fontsize=16)
    
    # Initialize each subplot
    lines = []
    for i in range(len(ylimits)):
        ax[i].set_ylim(ylimits[i])
        ax[i].set_xlim(xlimits)
        ax[i].set_ylabel(ylabels[i])
        ax[i].set_xlabel(xlabels[i])
        ax[i].grid(False)
        
        # Initialize lines within each subplot based on config
        line_handles = []
        for j in range(numMUs):
            line, = ax[i].plot([], [], color=colors[j % len(colors)], 
                               linestyle=lineStyles[j % len(lineStyles)], 
                               marker=markerStyles[j % len(markerStyles)], 
                               linewidth=linewidth, markersize=8)
            line_handles.append(line)
        lines.append(line_handles)

    plt.tight_layout()
    plt.show(block=False)
    
    return fig, ax, lines
