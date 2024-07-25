#   Authors: Christian Komo, Niels Bidault


import numpy as np
import matplotlib.pyplot as plt

#in plotDriftTransform, add legend and gausian distribution for x and y points
#Replace list variables that are unchanging with tuples, more efficient for calculations

#For functions like drifttransform, paramter should only be a single 2d array with all 6 initial variables
#Add legend for graphs 
#Same total pipe length but different interval should end with the same standard deviation


class beam:
    def __init__(self):
        pass

    # Ensure ellipse_polar has 'self' as the first parameter
    def ellipse_polar(self, t, a, b):
        return a * np.cos(t), b * np.sin(t)

    # Ensure methods that operate on instances have 'self' parameter
    def gen_6d_gaussian(self, mean, std_dev, num_particles=100):
        particles = np.random.normal(mean, std_dev, size=(num_particles, 6))
        return particles

    def std_6d(self, particles):
        std_devs = np.std(particles, axis=0)
        return tuple(std_devs)

    # Add 'self' to the method and use self.ellipse_polar to access the ellipse_polar method
    def plot_6d(self, values):
        
        position_x_values = values[:, 0]
        phase_x_values = values[:, 1]
        position_y_values = values[:, 2]
        phase_y_values = values[:, 3]
        energy_values = values[:, 4]
        time_values = values[:, 5]
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        std_dev_x = np.std(position_x_values)
        std_dev_phase_x = np.std(phase_x_values)
        std_dev_phase_y = np.std(phase_y_values)
        std_dev_y = np.std(position_y_values)
        std_dev_energy = np.std(energy_values)
        std_dev_time = np.std(time_values)

        plot_settings = [
            (position_x_values, phase_x_values, std_dev_x, std_dev_phase_x, 'Position x (mm)', 'Phase x (mrad)'),
            (position_y_values, phase_y_values, std_dev_y, std_dev_phase_y, 'Position y (mm)', 'Phase y (mrad)'),
            (position_x_values, position_y_values, std_dev_x, std_dev_y, 'Position x (mm)', 'Position y (mm)'),
            (energy_values, time_values, std_dev_energy, std_dev_time, 'Energy (keV)', 'Time (ns)'),
        ]

        for i, (x_vals, y_vals, a, b, xlabel, ylabel) in enumerate(plot_settings):
            ax = axes[i // 2, i % 2]

            t = np.linspace(0, 2 * np.pi, 100)
            # Use self.ellipse_polar to call the method
            ex, ey = self.ellipse_polar(t, a, b)

            ax.plot(ex, ey, color='black', linestyle='--', label='Ellipse')
            ax.plot(ex, -ey, color='black', linestyle='--')

            ax.scatter(x_vals, y_vals, s=15, alpha=0.7)

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.grid(True)

        plt.tight_layout()
        plt.show()
    
