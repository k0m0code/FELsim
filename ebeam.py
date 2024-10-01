#   Authors: Christian Komo, Niels Bidault


import numpy as np
import pandas as pd
import sympy as sp
from sympy.plotting import plot_implicit, PlotGrid
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from scipy.stats import norm


#in plotDriftTransform, add legend and gausian distribution for x and y points
#Replace list variables that are unchanging with tuples, more efficient for calculations
#Remove acceptance percentage when shapes are not in use
#Add legend for graphs 


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
        label_axes = ["x", "y", "z"]

        twiss_data = []

        for i in range(3):
            emittance = np.sqrt(
                dist_cov[2 * i, 2 * i] * dist_cov[2 * i + 1, 2 * i + 1] - dist_cov[2 * i, 2 * i + 1] ** 2)
            alpha = - dist_cov[2 * i, 2 * i + 1] / emittance
            beta = dist_cov[2 * i, 2 * i] / emittance
            gamma = dist_cov[2 * i + 1, 2 * i + 1] / emittance
            phi = 90 * np.arctan2(2 * alpha, gamma - beta) / np.pi

            twiss_data.append([emittance, alpha, beta, gamma, phi])

        twiss = pd.DataFrame(twiss_data, columns=label_twiss, index=label_axes[:3])

        return dist_avg, dist_cov, twiss

    def gen_6d_gaussian(self, mean, std_dev, num_particles=100):
        particles = np.random.normal(mean, std_dev, size=(num_particles, 6))
        return particles
    









    '''
    THIS FUNCTION BELOW IS A WIP
    '''

    def is_within_ellipse(self, x, y, xc, yc, twiss_axis, n):
        emittance = n * twiss_axis[r"$\epsilon$ ($\pi$.mm.mrad)"]
        alpha = twiss_axis[r"$\alpha$"]
        beta = twiss_axis[r"$\beta$ (m)"]
        gamma = twiss_axis[r"$\gamma$ (rad/m)"]

        # Calculate the ellipse equation
        Z = gamma * (x - xc) ** 2 + 2 * alpha * (x - xc) * (y - yc) + beta * (y - yc) ** 2 - emittance
        print("Z: ")
        print(Z)
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

          
            twiss_txt = '\n'.join(f'{label}: {np.round(value, 2)}' for label, value in twiss_axis.items())
            props = dict(boxstyle='round', facecolor='lavender', alpha=0.5)
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



















    


    '''
    returns twiss and eliptical data
    '''
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
        return std1, std6, dist_6d, twiss

    '''
    plots 6d and twiss data, used in schematic.py file
    '''
    def plotXYZ(self, dist_6d, std1, std6, twiss, ax1, ax2, ax3, ax4, maxVals = [0,0,0,0,0,0], minVals = [0,0,0,0,0,0], defineLim = True, shape = {}):
        axlist = [ax1,ax2,ax3]
        # Define SymPy symbols for plotting
        x_sym, y_sym = sp.symbols('x y')
        x_labels = [r'Position $x$ (mm)', r'Position $y$ (mm)', r'Time $\Delta$ $t$ (ns)', r'Position $x$ (mm)']
        y_labels = [r'Phase $x^{\prime}$ (mm)', r'Phase $y^{\prime}$ (mm)', r'Energy $\Delta$ $E$ (keV)',
                    r'Position $y$ (mm)']
        
        for i, axis in enumerate(['x', 'y', 'z']):
            twiss_axis = twiss.loc[axis]

            # Access Twiss parameters for the current axis
            ax = axlist[i]
            if defineLim:
                ax.set_xlim(minVals[2*i], maxVals[2*i])
                ax.set_ylim(minVals[2*i + 1], maxVals[2*i + 1])

            ax.scatter(dist_6d[:, 2 * i], dist_6d[:, 2 * i + 1], s=15, alpha=0.7)
            ax.contour(std1[i][0], std1[i][1], std1[i][2], levels=[0], colors='black', linestyles='--')
            ax.contour(std6[i][0], std6[i][1], std6[i][2], levels=[0], colors='black', linestyles='--')

            twiss_txt = '\n'.join(f'{label}: {np.round(value, 2)}' for label, value in twiss_axis.items())
            props = dict(boxstyle='round', facecolor='lavender', alpha=0.5)
            ax.text(0.01, 0.97, twiss_txt, transform=ax.transAxes, fontsize=8,
                    verticalalignment='top', bbox=props)

            ax.set_title(f'{axis} - Phase Space')
            ax.set_xlabel(x_labels[i])
            ax.set_ylabel(y_labels[i])
            ax.grid(True)

        # Plot for 'x, y - Space'
        withinArea = [[],[]]
        outsideArea = [[],[]]
        xyPart = [dist_6d[:, 0], dist_6d[:, 2]]
        if shape.get("shape") == "circle":
            radius = shape.get("radius")
            origin = shape.get("origin")
            for ii in range(len(xyPart[0])):
                x, y = xyPart[0][ii], xyPart[1][ii]
                if (x - origin[0])**2 + (y - origin[1])**2 < radius**2:
                    withinArea[0].append(x)
                    withinArea[1].append(y)
                else:
                    outsideArea[0].append(x)
                    outsideArea[1].append(y)
            circle = patches.Circle(origin, radius, edgecolor = "black", facecolor = "None")
            ax4.add_patch(circle)
        elif shape.get("shape") == "rectangle":
            length = shape.get("length")
            width = shape.get("width")
            origin = shape.get("origin")
            updatedOrigin = ((origin[0] - length/2),(origin[1] - width/2))
            for ii in range(len(xyPart[0])):
                x, y = xyPart[0][ii], xyPart[1][ii]
                if (origin[0] - (length/2)) < x < (origin[0] + (length/2)) and (origin[1] - (width/2)) < y < (origin[1] + (width/2)):
                    withinArea[0].append(x)
                    withinArea[1].append(y)
                else:
                    outsideArea[0].append(x)
                    outsideArea[1].append(y)
            rectangle = patches.Rectangle(updatedOrigin,length,width, edgecolor = "black", facecolor = "none")
            ax4.add_patch(rectangle)
        else:
            ax4.scatter(xyPart[0], xyPart[1], s=15,alpha = 0.7)
        if defineLim:
            ax4.set_xlim(minVals[0], maxVals[0])
            ax4.set_ylim(minVals[2], maxVals[2])
        ax4.scatter(withinArea[0], withinArea[1], s=15, alpha=0.7, color = "blue")
        ax4.scatter(outsideArea[0], outsideArea[1], s=15, alpha=0.7, color = "red")
        percentageInside = len(withinArea[0])/len(xyPart[0])*100
        

        ax4.set_title(f'x, y - Space: {round(percentageInside,4)}% acceptance')
        ax4.set_xlabel(x_labels[i + 1])
        ax4.set_ylabel(y_labels[i + 1])
        ax4.grid(True)

    '''
    plots 6d and twiss data with only particle distribution data
    '''
    def plot_6d(self, dist_6d, title):

        num_pts = 60  # Used for implicit plot of the ellipse
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

            # Plot the contour where Z = 0 (the ellipse)
            ax = axes[i // 2, i % 2]
            ax.scatter(dist_6d[:, 2 * i], dist_6d[:, 2 * i + 1], s=15, alpha=0.7)
            X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=1, num_pts=num_pts)
            ax.contour(X, Y, Z, levels=[0], colors='black', linestyles='--')
            X, Y, Z = self.ellipse_sym(dist_avg[2 * i], dist_avg[2 * i + 1], twiss_axis, n=6, num_pts=num_pts)
            ax.contour(X, Y, Z, levels=[0], colors='black', linestyles='--')

            # Construct the text string from the Twiss parameters
            twiss_txt = '\n'.join(f'{label}: {np.round(value, 2)}' for label, value in twiss_axis.items())
            props = dict(boxstyle='round', facecolor='lavender', alpha=0.5)
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
