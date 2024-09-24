#   Authors: Christian Komo, Niels Bidault

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
import csv
import numpy as np
from beamline import *
from ebeam import beam
import datetime
import time
from tqdm import tqdm

class draw_beamline:
    def __init__(self):
        '''
        Beamline object creates graphs and plots to visualize beamline data and information
        '''
        self.figsize = (10,9)  #  Size of graph window
        self.matrixVariables = None  # For access to data after transformation
        self.sixdValues = None  # For access to data after transformation


    '''
    POSSIBLE WIP FUNCTION IN THE FUTURE
    DISPLAY BEAMLINE AND DATA ABOUT EACH SEGMENT
    '''
    # def display_beamline(self, beamline):
    #     """
    #     Draws a schematic representation of a beamline.

    #     :param beamline: An object containing information about the beamline elements,
    #                      including their positions.
    #     """
    #     fig, ax = plt.subplots()

    #     # Set up the plot
    #     ax.set_xlim(min(beamline.z_start) - 1, max(beamline.z_end) + 1)
    #     ax.set_ylim(-1, 1)
    #     ax.set_xlabel('Position (m)')
    #     ax.set_title('Beamline Schematic')

    #     # Draw each beamline element
    #     for i, name in enumerate(beamline.names):
    #         z_start = beamline.z_start[i]
    #         z_end = beamline.z_end[i]
    #         element_width = z_end - z_start

    #         rect = patches.Rectangle((z_start, -self.element_height / 2), element_width, self.element_height,
    #                                  linewidth=self.stroke_width, edgecolor=self.element_color, facecolor='none')
    #         ax.add_patch(rect)
    #         ax.text((z_start + z_end) / 2, 0, name, ha='center', va='center')

    #     plt.grid(True)
    #     plt.show()



    def driftTransformScatter(self, values, length, plot = True):
        '''
        Simulates particles passing through drift space

        Parameters
        ----------
        values: np.array(list[float][float])
            2D numPy list of particle elements
        length: float
            length of the drift space particle passes through
        plot: bool, optional
            tells function whether to plot particle data or not

        Returns
        -------
        x_transform: list[float]
            list containing each particles' x position
        y_transform: list[float]
            list containing each particles' y position
        '''
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

    def csvWriteData(self, name, distance, std_x, std_y, mean_x, mean_y):
        '''
        Appends particle data to a csv file

        Parameters
        ----------
        name: str
            Name of csv file
        distance: float
            Length into the beamline element data is measured at
        std_x: float
            standard deviation measurement of x position
        std_y: float
            standard deviation measurement of y position
        mean_x: float
            average particle x position
        mean_y: float
            average particle y position
        '''
        with open(name, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            if csvfile.tell() == 0:
                csvwriter.writerow(['distance (mm)','x position standard deviaton (mm)', 'y position standard deviation (mm)', 'x position mean (mm)', 'y position mean (mm)'])
            csvwriter.writerow([distance, std_x, std_y, mean_x, mean_y])

    def checkMinMax(self, matrixVariables, maxval, minval):
        '''
        Updates max and min values of a set of particles in a beamline. Used for finding and
        setting the plot boundaries displaying the position of particles.

        Parameters
        ----------
        matrixvairables: np.array(list[float][float])
            A 6 column 2d numpy array, each row containing 6 initial values of each particle's measurements
        maxval: list[float]
            list of each current maximum value for each variable throughout the beamline
        minval: list[float]
            list of each current minimum value for each variable throughout the beamline

        Returns
        -------
        maxval: list[float]
            updated list of maximum values
        minval: list[float]
            updated list of minimum values
        '''
        initialx = matrixVariables[:, 0]
        initialy = matrixVariables[:, 2]
        initialxphase = matrixVariables[:, 1]
        initialyphase = matrixVariables[:, 3]
        initialz = matrixVariables[:, 4]
        initialzphase = matrixVariables[:, 5]
        initialList = [initialx, initialxphase, initialy, initialyphase, initialz, initialzphase]
        for i in range(len(initialList)):
            maximum = max(initialList[i])
            if maximum > maxval[i]:
                maxval[i] = maximum
            minimum = min(initialList[i])
            if minimum < minval[i]:
                minval[i] = minimum
        return maxval, minval

    def appendToList(self, xStd, yStd, xMean, yMean, x_axis, interval, matrixVariables):
        '''
        Append updated values to five different arrays, used for plotBeamPositionTransform

        Parameters
        ----------
        xStd: list[float]
            standard deviation of particles' x position for each distance interval
        yStd: list[float]
            standard deviation of particles' y position for each distance interval
        xMean: list[float]
            average of particles' x position for each distance interval
        yMean: list[float]
            average of particles' y position for each distance interval
        x_axis: list[float]
            contains distance intervals to measure particle data over
        interval: float
            the amount between each distance interval
        matrixVariables: np.array(list[float][float])
            2D numPy list of particle elements to measure data from

        Returns
        -------
        xStd: list[float]
            updated standard deviation of x position list
        yStd: list[float]
            updated standard deviation of y position list
        xMean: list[float]
            updated average of x position list
        yMean: list[float]
            updated average of y positiion list
        x[axis]: list[float]
            updated distance interval list
        '''
        xStd.append(np.std(matrixVariables[:,0]))
        yStd.append(np.std(matrixVariables[:,2]))
        xMean.append(np.mean(matrixVariables[:,0]))
        yMean.append(np.mean(matrixVariables[:,2]))
        x_axis.append(round(x_axis[-1]+interval, 3))
        return xStd, yStd, xMean, yMean, x_axis

    
    def plotBeamPositionTransform(self, matrixVariables, beamSegments, interval, defineLim = True, saveData = False, shape = {}, plot = True):
        '''
        Simulates movement of particles through an accelerator beamline

        Parameters
        ----------
        matrixvairables: np.array(list[float][float])
            A 6r x 10c 2d numpy array containing initial values of each electron's measurements
        beamSegmeents: list[beamline]
            Numpy array/list containing beamline objects which represent the beam
        interval: float
            Arbitrary number specifying interval for graph to take measurements at
        defineLim: bool, optional
            If plot should change dynamically with the points or stay static.
        saveData: boolean, optional
            Boolean value specifying whether to save data into a csv file or not
        shape: dict{}, optional
            dictionary storing info about the acceptance boundary
            ex. shape, width, radius, length, origin
        plot: bool, optional
            Optional boolean variable to plot simulation or not



        NOTE:
        shape is a dictionary defined as:
        shape = {"shape": "circle", "radius": 5, "origin": (0,5)}
        or
        shape = {"shape": "rectangle", "length": 200, "width": 500, "origin": (10,-4)}
        Only 2 shapes currently: rectangles and circles
        '''

        # Initialize values
        initialx = matrixVariables[:, 0]
        initialy = matrixVariables[:, 2]
        xStd = [np.std(initialx)]
        yStd = [np.std(initialy)]
        xMean = [np.mean(initialx)]
        yMean = [np.mean(initialy)]
        x_axis = [0]
        ebeam = beam()
        plot6dValues = {0: (ebeam.getXYZ(matrixVariables))}
        maxVals = [0, 0, 0, 0, 0, 0]
        minVals = [0, 0, 0, 0, 0, 0]

        if defineLim:
            maxVals, minVals = self.checkMinMax(matrixVariables, maxVals, minVals)

        total_intervals = sum(int(segment.length // interval) + 1 for segment in beamSegments)

        # Initialize the progress bar
        with tqdm(total=total_intervals, desc="Simulating Beamline",
                  bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]") as pbar:
            for i in range(len(beamSegments)):
                # Loop through each beamline object in beamSegments array
                intTrack = beamSegments[i].length

                while intTrack >= interval:
                    # Perform calculations to plot later on
                    # Use each segment's array to transform particles
                    matrixVariables = np.array(beamSegments[i].useMatrice(matrixVariables, length=interval))
                    xStd, yStd, xMean, yMean, x_axis = self.appendToList(xStd, yStd, xMean, yMean, x_axis, interval,
                                                                         matrixVariables)
                    if defineLim:
                        maxVals, minVals = self.checkMinMax(matrixVariables, maxVals, minVals)
                    plot6dValues.update({x_axis[-1]: (ebeam.getXYZ(matrixVariables))})

                    # Update the progress bar
                    pbar.update(1)
                    intTrack -= interval

                if intTrack > 0:
                    # Use remainder length once when it's smaller than interval
                    matrixVariables = np.array(beamSegments[i].useMatrice(matrixVariables, length=intTrack))
                    xStd, yStd, xMean, yMean, x_axis = self.appendToList(xStd, yStd, xMean, yMean, x_axis, intTrack,
                                                                         matrixVariables)
                    if defineLim:
                        maxVals, minVals = self.checkMinMax(matrixVariables, maxVals, minVals)
                    plot6dValues.update({x_axis[-1]: (ebeam.getXYZ(matrixVariables))})

                    # Update the progress bar for the remaining part
                    pbar.update(1)

        if saveData:
            #  Optionally save standard deviation and mean data
            name = "simulator-data-" + datetime.datetime.now().strftime('%Y-%m-%d') + "_" + datetime.datetime.now().strftime('%H_%M_%S') +".csv"
            for i in range(len(x_axis)):
                self.csvWriteData(name, x_axis[i], xStd[i], yStd[i], xMean[i], yMean[i])

        #  Testing purposes
        self.matrixVariables = matrixVariables
        self.sixdValues = plot6dValues

        if plot:
            #  Configure graph shape
            fig = plt.figure(figsize=self.figsize)
            gs = gridspec.GridSpec(3, 2, height_ratios=[0.8, 0.8, 1])
            ax1 = plt.subplot(gs[0,0])
            ax2 = plt.subplot(gs[0, 1])
            ax3 = plt.subplot(gs[1, 0])
            ax4 = plt.subplot(gs[1, 1])
            
            #  Plot inital 6d scatter data
            matrix = plot6dValues.get(0)
            ebeam.plotXYZ(matrix[2], matrix[0], matrix[1], matrix[3], ax1,ax2,ax3,ax4, maxVals, minVals, defineLim, shape)
        
            #  Plot and configure line graph data
            ax5 = plt.subplot(gs[2, :])
            plt.plot(x_axis, xStd, label = 'x position std')
            plt.plot(x_axis, yStd, label = "y position std")
            plt.plot(x_axis, xMean, color = 'red', label = 'x position mean')
            plt.plot(x_axis, yMean, color = 'blue', label = 'y position mean')
            ax5.set_xticks(x_axis)
            plt.xlim(0,x_axis[-1])
            ax5.set_xticklabels(x_axis,rotation=45,ha='right')
            plt.tick_params(labelsize = 9)
            plt.xlabel("Distance from start of beam (m)")
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
            
            #   Scroll bar creation and function
            scrollax = plt.axes([0.079,0.01,0.905,0.01], facecolor = 'lightgoldenrodyellow')
            scrollbar = Slider(scrollax, 'scroll', 0, x_axis[-1], valinit = 0, valstep=np.array(x_axis))
            def update_scroll(val):
                matrix = plot6dValues.get(scrollbar.val)
                ax1.clear()
                ax2.clear()
                ax3.clear()
                ax4.clear()
                ebeam.plotXYZ(matrix[2], matrix[0], matrix[1], matrix[3], ax1,ax2,ax3,ax4, maxVals, minVals, defineLim, shape)
                fig.canvas.draw_idle()
            scrollbar.on_changed(update_scroll)
            
            #  Title and ploting
            plt.suptitle("Beamline Simulation")
            plt.tight_layout()
            plt.show()
