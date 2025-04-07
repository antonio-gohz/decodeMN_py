import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg' depending on your system TkAgg
import matplotlib.pyplot as plt

# Default color palette
DEFAULT_COLORS = [
    (0.2667, 0.4667, 0.6667),
    (0.9333, 0.4000, 0.4667),
    (0.1333, 0.5333, 0.2000),
    (0.8, 0.7333, 0.2667),
    (0.4, 0.8, 0.9333),
    (0.6667, 0.2, 0.4667),
    (0.7333, 0.7333, 0.7333)
]

class LivePlot:
    def __init__(self, num_axes, num_lines=1):
        self.fig, self.ax = plt.subplots(num_axes, 1, figsize=(10, 7))
        self.lines = []
        self.num_axes = num_axes
        self.num_lines = num_lines
        self.is_open = True  # Track if the plot window is open

        
        # Initialize each subplot
        for i in range(num_axes):
            line_handles = []
            num_lines2 = num_lines[i] if i < len(num_lines) else 1  # Default to 1 if not specified
            for j in range(num_lines2):
                color = DEFAULT_COLORS[j % len(DEFAULT_COLORS)]
                line, = self.ax[i].plot([], [], color=color, linestyle='-', marker='none')
                line_handles.append(line)
            self.lines.append(line_handles)

        plt.tight_layout()
        plt.show(block=False)
        # Connect the figure close event to a handler
        self.fig.canvas.mpl_connect('close_event', self.on_close)

    def update(self, data):
        
        if not self.is_open:  # If the plot is closed, skip update
            return
        
        # data should be a list of lists: data[axis_index][line_index] = (x_data, y_data)
        for i, axis_data in enumerate(data):
            for j, (x_data, y_data) in enumerate(axis_data):
                self.lines[i][j].set_data(x_data, y_data)
                self.ax[i].relim()
                self.ax[i].autoscale_view()

        plt.draw()
        plt.pause(0.1)

    def on_close(self, event):
        """Handle the figure close event."""
        print("Figure closed!")
        self.is_open = False  # Mark the plot as closed when the figure is closed
        plt.close(self.fig)

    def close(self):
        """Gracefully close the plot."""
        plt.close(self.fig)
        self.is_open = False  # Mark the plot as closed
        print("Closing plot...")


# Example usage:
# Create a LivePlot instance with 2 axes and 1 line per axis
# live_plot = LivePlot(num_axes=2)

# To update the plots:
# live_plot.update([[(x_data1, y_data1)], [(x_data2, y_data2)]])
#     return x_data, y_data

# live_plot(data_generator)
