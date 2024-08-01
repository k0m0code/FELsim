#   Authors: Christian Komo, Niels Bidault


import numpy as np
import pandas as pd
import sympy as sp
from sympy.plotting import plot_implicit, PlotGrid
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#in plotDriftTransform, add legend and gausian distribution for x and y points
#Replace list variables that are unchanging with tuples, more efficient for calculations

#For functions like drifttransform, paramter should only be a single 2d array with all 6 initial variables
#Add legend for graphs 
#Same total pipe length but different interval should end with the same standard deviation


class beam:
    def __init__(self):
        pass

    x_sym, y_sym = sp.symbols('x_sym y_sym')

    def ellipse_sym(self, xc, yc, twiss_axis, n=1, num_pts=40):
        emittance = n * twiss_axis[r"$\epsilon$ ($\pi$.mm.mrad)"]
        alpha = twiss_axis[r"$\alpha$"]
        beta = twiss_axis[r"$\beta$ (m)"]
        gamma = twiss_axis[r"$\gamma$ (rad/m)"]

        x_max = xc + np.sqrt(emittance / (gamma - alpha ** 2 / beta))
        x_min = xc - np.sqrt(emittance / (gamma - alpha ** 2 / beta))
        y_max = yc + np.sqrt(emittance / (beta - alpha ** 2 / gamma))
        y_min = yc - np.sqrt(emittance / (beta - alpha ** 2 / gamma))

        x_vals = np.linspace(x_min, x_max, num_pts)
        y_vals = np.linspace(y_min, y_max, num_pts)
        X, Y = np.meshgrid(x_vals, y_vals)
        Z = gamma * (X - xc)** 2 + 2 * alpha * (X - xc) * (Y - yc) + beta * (Y - yc) ** 2 - emittance

        

        return X, Y, Z

    
    def cal_twiss(self, dist_6d, ddof=1):
        dist_avg = np.mean(dist_6d, axis=0)
        dist_cov = np.cov(dist_6d, rowvar=False, ddof=ddof)

        label_twiss = ["$\epsilon$ ($\pi$.mm.mrad)", r"$\alpha$", r"$\beta$ (m)", r"$\gamma$ (rad/m)", r"$\phi$ (deg)"]
        label_axes =["x", "y", "z"]
        twiss = pd.DataFrame(columns=label_twiss)
        for i in range(3):
            emittance = np.sqrt(dist_cov[2 * i, 2 * i] * dist_cov[2 * i + 1, 2 * i + 1] - dist_cov[2 * i, 2 * i + 1] ** 2)
            alpha = - dist_cov[2 * i, 2 * i + 1] / emittance
            beta = dist_cov[2 * i, 2 * i] / emittance
            gamma = dist_cov[2 * i + 1, 2 * i + 1] / emittance
            phi = 90 * np.arctan2(2 * alpha, gamma - beta) / np.pi
            tmp = pd.DataFrame([[emittance, alpha, beta, gamma, phi]], columns=label_twiss, index=[label_axes[i]])
            twiss = pd.concat([twiss, tmp])

        
        return dist_avg, dist_cov, twiss

    def gen_6d_gaussian(self, mean, std_dev, num_particles=100):
        particles = np.random.normal(mean, std_dev, size=(num_particles, 6))
        return particles
    









    '''
    THIS FUNCTION BELOW IS A WIP
    '''






    # def particles_in_ellipse(self, dist_6d, n = 1):
    #     fig, axes = plt.subplots(2, 2)

    #     dist_avg, dist_cov, twiss = self.cal_twiss(dist_6d, ddof=1)
    #     for i, axis in enumerate(['x', 'y', 'z']):
    #         ax = axes[i // 2, i % 2]
    #         twiss_axis = twiss.loc[axis]
    #         X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=n, num_pts=60)
    #         ax.scatter(dist_6d[:, 2 * i], dist_6d[:, 2 * i + 1], s=15, alpha=0.7)
    #         ax.contour(X, Y, Z, levels=[0], colors='black', linestyles='--')
    #         twiss_txt = '\n'.join(f'{label}: {np.round(value, 2)}' for label, value in twiss_axis.items())
    #         props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #         ax.text(0.05, 0.95, twiss_txt, transform=ax.transAxes, fontsize=12,
    #                 verticalalignment='top', bbox=props)

    #     plt.tight_layout()
    #     plt.show()
    #     return


    def is_within_ellipse(self, x, y, xc, yc, twiss_axis, n):
        emittance = n * twiss_axis[r"$\epsilon$ ($\pi$.mm.mrad)"]
        alpha = twiss_axis[r"$\alpha$"]
        beta = twiss_axis[r"$\beta$ (m)"]
        gamma = twiss_axis[r"$\gamma$ (rad/m)"]

        # Calculate the ellipse equation
        Z = gamma * (x - xc) ** 2 + 2 * alpha * (x - xc) * (y - yc) + beta * (y - yc) ** 2 - emittance
        # print("Z: ")
        # print(Z)
        # Check if the point (x, y) is inside or on the ellipse
        return Z <= 0


    def particles_in_ellipse(self, dist_6d, n=1):
        fig, axes = plt.subplots(2, 2)

        dist_avg, dist_cov, twiss = self.cal_twiss(dist_6d, ddof=1)
        count_within_ellipse = []  # List to store count of particles within ellipse for each axis

        for i, axis in enumerate(['x', 'y', 'z']):
            ax = axes[i // 2, i % 2]
            twiss_axis = twiss.loc[axis]
            X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=n, num_pts=60)

            # Plot particle points and ellipse contours
            ax.scatter(dist_6d[:, 2 * i], dist_6d[:, 2 * i + 1], s=15, alpha=0.7)
            ax.contour(X, Y, Z, levels=[0], colors='black', linestyles='--')

            # Evaluate particles within ellipse
            num_within_ellipse = 0 
            for j in range(len(dist_6d)):
                x, y = dist_6d[j, 2 * i], dist_6d[j, 2 * i + 1] #  Retrives x and y from each value list
                # Check if the particle is within the ellipse
                if self.is_within_ellipse(x, y, dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n = n):
                    num_within_ellipse += 1

            # num_within_ellipse = 0
            # for j in Z:
            #     for jj in j:
            #         if jj <= 0:
            #             num_within_ellipse += 1

            twiss_txt = '\n'.join(f'{label}: {np.round(value, 2)}' for label, value in twiss_axis.items())
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.05, 0.95, twiss_txt, transform=ax.transAxes, fontsize=12,
                    verticalalignment='top', bbox=props)
            
            count_within_ellipse.append(num_within_ellipse)
            # Optional: Display the count on the plot
            ax.set_title(f'{axis} - Phase Space\nParticles within ellipse: {num_within_ellipse}')
            ax.set_xlabel(f'Position {axis}')
            ax.set_ylabel(f'Phase {axis} prime')
            ax.grid(True)

        plt.tight_layout()
        plt.show()

        return count_within_ellipse






















    # def plot_6d(self, dist_6d, title):

    #     num_pts = 60  # Used for implicit plot of the ellipse
    #     ddof = 1  # Unbiased Bessel correction for standard deviation calculation

    #     dist_avg, dist_cov, twiss = self.cal_twiss(dist_6d, ddof=ddof)

    #     # Define SymPy symbols for plotting
    #     x_sym, y_sym = sp.symbols('x y')

    #     fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
    #     x_labels = [r'Position $x$ (mm)', r'Position $y$ (mm)', r'Energy $\Delta$ $E$ (keV)', r'Position $x$ (mm)']
    #     y_labels = [r'Phase $x^{\prime}$ (mm)', r'Phase $y^{\prime}$ (mm)', r'Time $\Delta$ $t$ (ns)',
    #                 r'Position $y$ (mm)']

    #     for i, axis in enumerate(['x', 'y', 'z']):

    #         # Access Twiss parameters for the current axis
    #         twiss_axis = twiss.loc[axis]

    #         # Plot the contour where Z = 0 (the ellipse)
    #         ax = axes[i // 2, i % 2]
    #         ax.scatter(dist_6d[:, 2 * i], dist_6d[:, 2 * i + 1], s=15, alpha=0.7)
    #         X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=1, num_pts=num_pts)
    #         ax.contour(X, Y, Z, levels=[0], colors='black', linestyles='--')
    #         X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=6, num_pts=num_pts)
    #         ax.contour(X, Y, Z, levels=[0], colors='black', linestyles='--')

    #         # Construct the text string from the Twiss parameters
    #         twiss_txt = '\n'.join(f'{label}: {np.round(value, 2)}' for label, value in twiss_axis.items())
    #         props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #         ax.text(0.05, 0.95, twiss_txt, transform=ax.transAxes, fontsize=12,
    #                 verticalalignment='top', bbox=props)

    #         ax.set_title(f'{axis} - Phase Space')
    #         ax.set_xlabel(x_labels[i])
    #         ax.set_ylabel(y_labels[i])
    #         ax.grid(True)

    #     # Plot for 'x, y - Space'
    #     ax = axes[(i + 1) // 2, (i + 1) % 2]
    #     ax.scatter(dist_6d[:, 0], dist_6d[:, 2], s=15, alpha=0.7)

    #     ax.set_title(f'x, y - Space')
    #     ax.set_xlabel(x_labels[i + 1])
    #     ax.set_ylabel(y_labels[i + 1])
    #     ax.grid(True)

    #     plt.suptitle(title)
    #     plt.tight_layout()
    #     plt.show()

    def getXYZ(self, dist_6d):
        num_pts = 60  # Used for implicit plot of the ellipse
        ddof = 1  # Unbiased Bessel correction for standard deviation calculation
        dist_avg, dist_cov, twiss = self.cal_twiss(dist_6d, ddof=ddof)
        std6 = []
        std1 = []
        for i, axis in enumerate(['x', 'y', 'z']):
            # Access Twiss parameters for the current axis
            twiss_axis = twiss.loc[axis]
            X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=1, num_pts=num_pts)
            std1.append([X,Y,Z])
            X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=6, num_pts=num_pts)
            std6.append([X,Y,Z])
        return std1, std6
    

    def plotXYZ(self, dist_6d, std1, std6, ax1, ax2, ax3, ax4):
        num_pts = 60  # Used for implicit plot of the ellipse
        axlist = [ax1,ax2,ax3]
        # Define SymPy symbols for plotting
        x_sym, y_sym = sp.symbols('x y')
        x_labels = [r'Position $x$ (mm)', r'Position $y$ (mm)', r'Energy $\Delta$ $E$ (keV)', r'Position $x$ (mm)']
        y_labels = [r'Phase $x^{\prime}$ (mm)', r'Phase $y^{\prime}$ (mm)', r'Time $\Delta$ $t$ (ns)',
                    r'Position $y$ (mm)']
        dist_avg, dist_cov, twiss = self.cal_twiss(dist_6d, ddof=ddof)
        
        for i, axis in enumerate(['x', 'y', 'z']):
            twiss_axis = twiss.loc[axis]

            # Access Twiss parameters for the current axis
            ax = axlist[i]
            ddof = 1  # Unbiased Bessel correction for standard deviation calculation

            ax.scatter(dist_6d[:, 2 * i], dist_6d[:, 2 * i + 1], s=15, alpha=0.7)
            ax.contour(std1[0], std1[1], std1[2], levels=[0], colors='black', linestyles='--')
            ax.contour(std6[0], std6[1], std6[2], levels=[0], colors='black', linestyles='--')

            twiss_txt = '\n'.join(f'{label}: {np.round(value, 2)}' for label, value in twiss_axis.items())
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.05, 0.95, twiss_txt, transform=ax.transAxes, fontsize=12,
                    verticalalignment='top', bbox=props)

            ax.set_title(f'{axis} - Phase Space')
            ax.set_xlabel(x_labels[i])
            ax.set_ylabel(y_labels[i])
            ax.grid(True)

        # Plot for 'x, y - Space'
        ax4.scatter(dist_6d[:, 0], dist_6d[:, 2], s=15, alpha=0.7)

        ax4.set_title(f'x, y - Space')
        ax4.set_xlabel(x_labels[i + 1])
        ax4.set_ylabel(y_labels[i + 1])
        ax4.grid(True)



    def plot_6d(self, dist_6d, ax1, ax2, ax3, ax4):

        num_pts = 60  # Used for implicit plot of the ellipse
        ddof = 1  # Unbiased Bessel correction for standard deviation calculation

        dist_avg, dist_cov, twiss = self.cal_twiss(dist_6d, ddof=ddof)
        axlist = [ax1,ax2,ax3]

        # Define SymPy symbols for plotting
        x_sym, y_sym = sp.symbols('x y')
        
        x_labels = [r'Position $x$ (mm)', r'Position $y$ (mm)', r'Energy $\Delta$ $E$ (keV)', r'Position $x$ (mm)']
        y_labels = [r'Phase $x^{\prime}$ (mm)', r'Phase $y^{\prime}$ (mm)', r'Time $\Delta$ $t$ (ns)',
                    r'Position $y$ (mm)']

        for i, axis in enumerate(['x', 'y', 'z']):

            # Access Twiss parameters for the current axis
            twiss_axis = twiss.loc[axis]

            # Plot the contour where Z = 0 (the ellipse)
            ax = axlist[i]
            ax.scatter(dist_6d[:, 2 * i], dist_6d[:, 2 * i + 1], s=15, alpha=0.7)
            X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=1, num_pts=num_pts)
            ax.contour(X, Y, Z, levels=[0], colors='black', linestyles='--')
            X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=6, num_pts=num_pts)
            ax.contour(X, Y, Z, levels=[0], colors='black', linestyles='--')

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
        ax4.scatter(dist_6d[:, 0], dist_6d[:, 2], s=15, alpha=0.7)

        ax4.set_title(f'x, y - Space')
        ax4.set_xlabel(x_labels[i + 1])
        ax4.set_ylabel(y_labels[i + 1])
        ax4.grid(True)
