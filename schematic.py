# schematic.py

import matplotlib.pyplot as plt
import matplotlib.patches as patches

class draw_beamline:
    def __init__(self):
        self.stroke_width = 2  # Example of a constant related to the figure

        # Other figure-related constants can be defined here
        self.element_height = 0.2
        self.element_color = 'blue'

    def display_beamline(self, beamline):
        """
        Draws a schematic representation of a beamline.

        :param beamline: An object containing information about the beamline elements,
                         including their positions.
        """
        fig, ax = plt.subplots()

        # Set up the plot
        ax.set_xlim(min(beamline.z_start) - 1, max(beamline.z_end) + 1)
        ax.set_ylim(-1, 1)
        ax.set_xlabel('Position (m)')
        ax.set_title('Beamline Schematic')

        # Draw each beamline element
        for i, name in enumerate(beamline.names):
            z_start = beamline.z_start[i]
            z_end = beamline.z_end[i]
            element_width = z_end - z_start

            rect = patches.Rectangle((z_start, -self.element_height / 2), element_width, self.element_height,
                                     linewidth=self.stroke_width, edgecolor=self.element_color, facecolor='none')
            ax.add_patch(rect)
            ax.text((z_start + z_end) / 2, 0, name, ha='center', va='center')

        plt.grid(True)
        plt.show()
