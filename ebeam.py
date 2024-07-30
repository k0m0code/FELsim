#   Authors: Christian Komo, Niels Bidault


import numpy as np
import pandas as pd
import sympy as sp
from sympy.plotting import plot_implicit, PlotGrid
import matplotlib.pyplot as plt


#in plotDriftTransform, add legend and gausian distribution for x and y points
#Replace list variables that are unchanging with tuples, more efficient for calculations

#For functions like drifttransform, paramter should only be a single 2d array with all 6 initial variables
#Add legend for graphs 
#Same total pipe length but different interval should end with the same standard deviation


class beam:
    def __init__(self):
        pass

    x_sym, y_sym = sp.symbols('x_sym y_sym')

    def ellipse_sym(self, x_sym, y_sym, xc, yc, twiss_axis):
        emittance = twiss_axis[r"$\epsilon$ ($\pi$.mm.mrad)"]
        alpha = twiss_axis[r"$\alpha$"]
        beta = twiss_axis[r"$\beta$ (m)"]
        gamma = twiss_axis[r"$\gamma$ (rad/m)"]
        return gamma * (x_sym - xc)** 2 + 2 * alpha * (x_sym - xc) * (y_sym - yc) + beta * (y_sym - yc) ** 2 - emittance

    def particles_in_ellipse(self):
        return

    def cal_twiss(self, dist_6d, ddof=1):
        dist_avg = np.mean(dist_6d, axis=0)
        dist_cov = np.cov(dist_6d, rowvar=False, ddof=ddof)

        label_twiss = ["$\epsilon$ ($\pi$.mm.mrad)", r"$\alpha$", r"$\beta$ (m)", r"$\gamma$ (rad/m)", "Phi (deg)"]
        label_axes =["x", "y", "z"]
        twiss = pd.DataFrame(columns=label_twiss)
        for i in range(3):
            emittance = np.sqrt(dist_cov[2 * i, 2 * i] * dist_cov[2 * i + 1, 2 * i + 1] - dist_cov[2 * i, 2 * i + 1] ** 2)
            alpha = - dist_cov[2 * i, 2 * i + 1] / emittance
            beta = dist_cov[2 * i, 2 * i] / emittance
            gamma = dist_cov[2 * i + 1, 2 * i + 1] / emittance
            phi = 180 / np.pi * np.arctan2(2 * alpha, gamma - beta) / 2
            tmp = pd.DataFrame([[emittance, alpha, beta, gamma, phi]], columns=label_twiss, index=[label_axes[i]])
            twiss = pd.concat([twiss, tmp])
        return dist_avg, dist_cov, twiss

    def gen_6d_gaussian(self, mean, std_dev, num_particles=100):
        particles = np.random.normal(mean, std_dev, size=(num_particles, 6))
        return particles

    def plot_6d(self, dist_6d, title):
        import matplotlib.pyplot as plt
        import sympy as sp
        import numpy as np

        num_pts = 40  # Used for implicit plot of the ellipse
        ddof = 1  # Unbiased Bessel correction for standard deviation calculation

        dist_avg, dist_cov, twiss = self.cal_twiss(dist_6d, ddof=ddof)

        # Define SymPy symbols for plotting
        x_sym, y_sym = sp.symbols('x y')

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        x_labels = [r'Position $x$ (mm)', r'Position $y$ (mm)', r'Energy $\Delta$ $E$ (keV)', r'Position $x$ (mm)']
        y_labels = [r'Phase $x^{\prime}$ (mm)', r'Phase $y^{\prime}$ (mm)', r'Time $\Delta$ $t$ (ns)',
                    r'Position $y$ (mm)']

        for i, axis in enumerate(['x', 'y', 'z']):
            # Access Twiss parameters for the current axis
            twiss_axis = twiss.loc[axis]

            # Define the ellipse equation using SymPy
            ellipse_eq = self.ellipse_sym(x_sym, y_sym, dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis)

            # Create a lambda function from the SymPy expression
            ellipse_func = sp.lambdify((x_sym, y_sym), ellipse_eq, modules='numpy')

            # Generate a grid of points
            emittance = twiss_axis[r"$\epsilon$ ($\pi$.mm.mrad)"]
            alpha = twiss_axis[r"$\alpha$"]
            beta = twiss_axis[r"$\beta$ (m)"]
            gamma = twiss_axis[r"$\gamma$ (rad/m)"]

            x_max = dist_avg[2 * i] + np.sqrt(emittance / (gamma - alpha ** 2 / beta))
            x_min = dist_avg[2 * i] - np.sqrt(emittance / (gamma - alpha ** 2 / beta))
            y_max = dist_avg[2 * i + 1] + np.sqrt(emittance / (beta - alpha ** 2 / gamma))
            y_min = dist_avg[2 * i + 1] - np.sqrt(emittance / (beta - alpha ** 2 / gamma))

            x_vals = np.linspace(x_min, x_max, num_pts)
            y_vals = np.linspace(y_min, y_max, num_pts)
            X, Y = np.meshgrid(x_vals, y_vals)
            Z = ellipse_func(X, Y)

            # Plot the contour where Z = 0 (the ellipse)
            ax = axes[i // 2, i % 2]
            ax.scatter(dist_6d[:, 2 * i], dist_6d[:, 2 * i + 1], s=15, alpha=0.7)
            ax.contour(X, Y, Z, levels=[0], colors='black')

            # Construct the text string from the Twiss parameters
            twiss_txt = '\n'.join(f'{label}: {np.round(value, 2)}' for label, value in twiss_axis.items())
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.05, 0.95, twiss_txt, transform=ax.transAxes, fontsize=12,
                    verticalalignment='top', bbox=props)

            ax.set_title(f'{axis} - Phase Space')
            ax.set_xlabel(x_labels[i])
            ax.set_ylabel(y_labels[i])
            ax.grid(True)

        # Plot for 'x, y - Space'
        ax = axes[(i + 1) // 2, (i + 1) % 2]
        ax.scatter(dist_6d[:, 0], dist_6d[:, 2], s=15, alpha=0.7)

        ax.set_title(f'x, y - Space')
        ax.set_xlabel(x_labels[i + 1])
        ax.set_ylabel(y_labels[i + 1])
        ax.grid(True)

        plt.suptitle(title)
        plt.tight_layout()
        plt.show()
