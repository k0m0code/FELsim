#   Authors: Christian Komo, Niels Bidault

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
from scipy.stats import norm
import csv
import numpy as np
from beamline import *
from ebeam import beam
import datetime

#Same total pipe length but different interval should end with the same standard deviation (PROBLEMS WITH QPF AND QPD FUNC)
#Make spacing of chart more efficient and useful?

class draw_beamline:
    def __init__(self):
        self.stroke_width = 2  # Example of a constant related to the figure

        # Other figure-related constants can be defined here
        self.element_height = 0.2
        self.element_color = 'blue'


    def csvWriteData(self, name, distance, std_x, std_y, mean_x, mean_y):
        with open(name, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            if csvfile.tell() == 0:
                csvwriter.writerow(['distance (mm)','x position standard deviaton (mm)', 'y position standard deviation (mm)', 'x position mean (mm)', 'y position mean (mm)'])
            csvwriter.writerow([distance, std_x, std_y, mean_x, mean_y])


    '''
    Class function to append updated values to five different arrays
    '''
    def appendToList(self, xStd, yStd, xMean, yMean, x_axis, interval, matrixVariables):
        xStd.append(np.std(matrixVariables[:,0]))
        yStd.append(np.std(matrixVariables[:,2]))
        xMean.append(np.mean(matrixVariables[:,0]))
        yMean.append(np.mean(matrixVariables[:,2]))
        x_axis.append(round(x_axis[-1]+interval, 3))
        return xStd, yStd, xMean, yMean, x_axis

    '''
    matrixvairables: np.array(list[float][float])
    A 6r x 10c 2d numpy array containing initial values of each electron's measurements

    beamSegmeents: list[beamline]
    numpy array/list containing beamline objects which represent the beam

    interval: float
    arbitrary number specifying interval for graph to take measurements at

    plot_z: tuple(float)
    tuple of numbered points in the beam on where to plot 6D

    saveData: boolean
    boolean value specifying whether to save data into a csv file or not
    '''
    def plotBeamPositionTransform(self, matrixVariables, beamSegments, interval, saveData = False):
        #  Initialize values
        xStd = [np.std(matrixVariables[:, 0])]
        yStd = [np.std(matrixVariables[:, 2])]
        xMean = [np.mean(matrixVariables[:, 0])]
        yMean = [np.mean(matrixVariables[:, 2])]
        x_axis = [0]
        xaxisMax = 0
        plot6dValues = {0: matrixVariables}

        #  Loop through each beamline object in beamSegments array
        for i in range(len(beamSegments)):
            intTrack = beamSegments[i].length
            xaxisMax += intTrack

            #  Make calculations to plot later on
            while intTrack >= interval:
                matrixVariables = np.array(beamSegments[i].useMatrice(matrixVariables, length = interval))
                xStd, yStd, xMean, yMean, x_axis = self.appendToList(xStd, yStd, xMean, yMean, x_axis, interval, matrixVariables)
                plot6dValues.update({x_axis[-1]: matrixVariables})
                intTrack -= interval
            if intTrack > 0:
                matrixVariables = np.array(beamSegments[i].useMatrice(matrixVariables, length = intTrack))
                xStd, yStd, xMean, yMean, x_axis = self.appendToList(xStd, yStd, xMean, yMean, x_axis, intTrack, matrixVariables)
                plot6dValues.update({x_axis[-1]: matrixVariables})

        #  Optionally save standard deviation and mean data
        if saveData:
            name = "simulator-data-" + datetime.datetime.now().strftime('%Y-%m-%d') + "_" + datetime.datetime.now().strftime('%H_%M_%S') +".csv"
            for i in range(len(x_axis)):
                self.csvWriteData(name, x_axis[i], xStd[i], yStd[i], xMean[i], yMean[i])

        #  Configure graph shape
        fig = plt.figure(figsize=(10, 9))
        gs = gridspec.GridSpec(3, 2, height_ratios=[0.8, 0.8, 1])
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0, 1])
        ax3 = plt.subplot(gs[1, 0])
        ax4 = plt.subplot(gs[1, 1])
        
        #  Plot inital 6d scatter data
        ebeam = beam()
        ebeam.plot_6d(matrixVariables, ax1,ax2,ax3,ax4)
        
        #  Plot and optimize line graph data
        ax5 = plt.subplot(gs[2, :])
        plt.plot(x_axis, xStd, label = 'x position std')
        plt.plot(x_axis, yStd, label = "y position std")
        plt.plot(x_axis, xMean, color = 'red', label = 'x position mean')
        plt.plot(x_axis, yMean, color = 'blue', label = 'y position mean')
        ax5.set_xticks(x_axis)
        plt.xlim(0,x_axis[-1])
        ax5.set_xticklabels(x_axis,rotation=45,ha='right')
        plt.tick_params(labelsize = 9)
        plt.xlabel("Distance from start of beam (mm)")
        plt.ylabel("Standard deviation (mm)")
        plt.legend()

        #  Create visual representation of beamline segments
        ymin, ymax = ax5.get_ylim()
        ax5.set_ylim(ymin-(ymax*0.05), ymax)
        ymin, ymax = ax5.get_ylim()
        blockstart = 0
        for seg in beamSegments:
            rectangle = patches.Rectangle((blockstart, ymin), seg.length, ymax*0.05, linewidth=1, edgecolor=seg.color, facecolor= seg.color)
            ax5.add_patch(rectangle)
            blockstart += seg.length
        
        scrollax = plt.axes([0.079,0.01,0.905,0.01], facecolor = 'lightgoldenrodyellow')
        scrollbar = Slider(scrollax, 'scroll', 0, x_axis[-1], valinit = 0, valstep=np.array(x_axis))
        def update_scroll(val):
            matrix = plot6dValues.get(scrollbar.val)
            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax4.clear()
            ebeam.plot_6d(matrix, ax1,ax2,ax3,ax4)
            fig.canvas.draw_idle()
        scrollbar.on_changed(update_scroll)
        
        plt.suptitle("Beamline Simulation")
        plt.tight_layout()
        plt.show()


    



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



    # Can length variable be negative?
    def driftTransformScatter(self, values, length, plot = True):

        x_pos = values[:, 0]
        phase_x = values[:, 1]
        y_pos = values[:, 2]
        phase_y = values[:, 3]
        x_transform = []
        y_transform = []

        for i in range(len(x_pos)):
            x_transform.append(x_pos[i]+length*phase_x[i])
        for i in range(len(y_pos)):
            y_transform.append(y_pos[i] + length*phase_y[i])

        if plot:
            fig, ax = plt.subplots()
            ax.scatter(x_pos,y_pos, c = 'blue', s=15, alpha=0.7, label = "Initial values")
            ax.scatter(x_transform, y_transform, c = 'green', s=15, alpha=0.7, label = "Transformed values")
            ax.set_xlabel('Position x (mm)')
            ax.set_ylabel('Position y (mm)')
            plt.legend(loc = 'upper right')
            plt.suptitle("Drift Transformation over " + str(length) + " mm")
            plt.tight_layout()
            plt.show()

        return x_transform, y_transform
