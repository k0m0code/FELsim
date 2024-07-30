#   Authors: Christian Komo, Niels Bidault

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import csv
import numpy as np
from beamline import *
from ebeam import beam
import datetime
from matplotlib.widgets import Slider

#Add legend for graphs 
#Same total pipe length but different interval should end with the same standard deviation

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

    plot_z: float
    arbitrary number specifying point in beamline to plot values in 6D

    saveData: boolean
    boolean value specifying whether to save data into a csv file or not
    '''
    def plotBeamPositionTransform(self, matrixVariables, beamSegments, interval, plot_z = -1, saveData = False):
        #  Initialize values
        xStd = [np.std(matrixVariables[:, 0])]
        yStd = [np.std(matrixVariables[:, 2])]
        xMean = [np.mean(matrixVariables[:, 0])]
        yMean = [np.mean(matrixVariables[:, 2])]
        x_axis = [0]
        xaxisMax = 0

        #  Loop through each beamline object in beamSegments array
        for i in range(len(beamSegments)):
            intTrack = beamSegments[i].length
            xaxisMax += intTrack

            #  Check if our z point is within the interval of our beamline object
            if ((xaxisMax - intTrack) < plot_z <= xaxisMax):
                sixDarry = np.array(beamSegments[i].useMatrice(matrixVariables, length = (plot_z-(xaxisMax-intTrack))))
                ebeam = beam()
                ebeam.plot_6d(sixDarry)

            #  Make calculations to plot later on
            while intTrack >= interval:
                matrixVariables = np.array(beamSegments[i].useMatrice(matrixVariables, length = interval))
                xStd, yStd, xMean, yMean, x_axis = self.appendToList(xStd, yStd, xMean, yMean, x_axis, interval, matrixVariables)
                intTrack -= interval
            if intTrack > 0:
                matrixVariables = np.array(beamSegments[i].useMatrice(matrixVariables, length = intTrack))
                xStd, yStd, xMean, yMean, x_axis = self.appendToList(xStd, yStd, xMean, yMean, x_axis, intTrack, matrixVariables)

        if saveData:
            name = "simulator-data-" + datetime.datetime.now().strftime('%Y-%m-%d') + "_" + datetime.datetime.now().strftime('%H_%M_%S') +".csv"
            for i in range(len(x_axis)):
                self.csvWriteData(name, x_axis[i], xStd[i], yStd[i], xMean[i], yMean[i])

        # print("x: \n" + str(x_axis))
        # print("xstd: \n" + str(xStd))
        # print ("ystd: \n" + str(yStd))

        fig, ax = plt.subplots()
        plt.plot(x_axis, xStd, label = 'x position std')
        plt.plot(x_axis, yStd, label = "y position std")
        plt.plot(x_axis, xMean, color = 'red', label = 'x position mean')
        plt.plot(x_axis, yMean, color = 'blue', label = 'y position mean')

        plt.subplots_adjust(bottom=0.25)
        ax.set_xticks(x_axis)
        plt.xlim(0,x_axis[-1] + x_axis[-1]*0.10)
        ax.set_xticklabels(x_axis,rotation=45,ha='right')
        plt.tick_params(labelsize = 9)
        plt.xlabel("Distance from start of beam (mm)")
        plt.ylabel("Standard deviation (mm)")
        plt.xlim(0, xaxisMax*0.20)
        plt.legend()

        scrollax = plt.axes([0.1,0.02,0.8,0.06], facecolor = 'lightgoldenrodyellow')
        scrollbar = Slider(scrollax, 'scroll', 0, 100, valinit = 0, valstep=1)
        def update_scroll(val):
            pos = scrollbar.val
            ax.set_xlim((pos/125)*xaxisMax, (pos/125)*xaxisMax + xaxisMax*(0.2))
            fig.canvas.draw_idle()
        scrollbar.on_changed(update_scroll)
        
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



    #integrate getDriftMatrice within this function or get rid of getDriftMatrice altogether
    # Can length variable be negative?
    def driftTransformScatter(self, x_pos: list[int], y_pos: list[int], phase_x: list[int], phase_y: list[int], length, plot = True):
        x_transform = []
        y_transform = []

        for i in range(len(x_pos)):
            x_transform.append(x_pos[i]+length*phase_x[i])
        for i in range(len(y_pos)):
            y_transform.append(y_pos[i] + length*phase_y[i])

        if plot:
            fig, ax = plt.subplots()
            ax.scatter(x_pos,y_pos, c = 'blue', s=15, alpha=0.7)
            ax.scatter(x_transform, y_transform, c = 'green', s=15, alpha=0.7)
            ax.set_xlabel('Position x (mm)')
            ax.set_ylabel('Position y (mm)')
            plt.legend(loc = 'upper right')
            plt.tight_layout()
            plt.show()

        return x_transform, y_transform
