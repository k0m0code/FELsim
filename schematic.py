#   Authors: Christian Komo, Niels Bidault

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
import pandas as pd
import csv
import numpy as np
from beamline import *
from ebeam import beam
import datetime
from tqdm import tqdm
from matplotlib.widgets import Button

#Note: plotBeamTransform may show rounded interval values slightly innacurate to actual lengths, but calculations are made with exact values, (rounded values only used for visualization)

class draw_beamline:
    def __init__(self):
        '''
        Beamline object creates graphs and plots to visualize beamline data and information
        '''
        self.figsize = (10,9)  #  Size of graph window
        self.matrixVariables = None  # For access to data after transformation
        self.sixdValues = None  # For access to data after transformation
        self.DEFAULTINTERVAL = 0.05
        self.DEFAULTINTERVALROUND = 5  # Number of decimal places p
        self.DEFAULTSPACINGPERCENTAGE = 0.02  # percent of total beamline length before adding plot tick


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

    

    def _createLabels(self, xaxis, spacing):
        xLabels = []

        defaultSpace = spacing
        for i in range(len(xaxis)):
            if i % defaultSpace == 0:
                xLabels.append(xaxis[i])
            else:
                xLabels.append("")
        return xLabels


    
    def _csvWriteData(self, name, twiss, x_axis):
        with open(name, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            if csvfile.tell() == 0:
                labels = ["z distance"]
                for lab in twiss:
                    for axis in twiss.index:
                        labels.append(str(axis) + ": " + str(lab))
                csvwriter.writerow(labels)
            for i in range(len(x_axis)):
                data = [x_axis[i]]
                for lab in twiss:
                        for axis in twiss.index:
                            list = twiss.at[axis, lab]
                            data.append(list[i])
                csvwriter.writerow(data)
                      
    def setEqualAxisScaling(self, maxVals, minVals):
        if maxVals[0] > maxVals[2]:
            maxVals[2] = maxVals[0]
        else:
            maxVals[0] = maxVals[2]
        if maxVals[1] > maxVals[3]:
            maxVals[3] = maxVals[1]
        else:
            maxVals[1] = maxVals[3]
        if minVals[0] < minVals[2]:
            minVals[2] = minVals[0]
        else:
            minVals[0] = minVals[2]
        if minVals[1] < minVals[3]:
            minVals[3] = minVals[1]
        else:
            minVals[1] = minVals[3]

    def plotBeamPositionTransform(self, matrixVariables, beamSegments, interval = -1, defineLim = True,
                                   saveData = False, shape = {}, plot = True, spacing = True, matchScaling = True):
        '''
        Simulates movement of particles through an accelerator beamline

        Parameters
        ----------
        matrixvairables: np.array(list[float][float])
            A 6r x 10c 2d numpy array containing initial values of each electron's measurements
        beamSegments: list[beamline]
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
        spacing: bool, optional
            Optional variable to optimize spacing of x labels when plotting for readibility
        matchScaling: bool, optional
            Whether to have same x' vs x and y' vs y axis scaling or not.
            defineLim must be True for same scaling setting to work



        NOTE:
        shape is a dictionary defined as:
        shape = {"shape": "circle", "radius": 5, "origin": (0,5)}
        or
        shape = {"shape": "rectangle", "length": 200, "width": 500, "origin": (10,-4)}
        Only 2 shapes currently: rectangles and circles
        '''

        # Initialize values
        ebeam = beam()
        result = ebeam.getXYZ(matrixVariables)
        twiss = result[3]
        plot6dValues = {0: result}
        twiss_aggregated_df = pd.DataFrame(
            {axis: {label: [] for label in twiss.index} for axis in twiss.columns}
        )  # Create twiss dataframe for particles
        # Add initial twiss values
        for i, axis in enumerate(twiss.index):
            twiss_axis = twiss.loc[axis]
            for label, value in twiss_axis.items():
                twiss_aggregated_df.at[axis, label].append(value)

        x_axis = [0]
        maxVals = [0, 0, 0, 0, 0, 0]
        minVals = [0, 0, 0, 0, 0, 0]

        if defineLim:
            maxVals, minVals = self.checkMinMax(matrixVariables, maxVals, minVals)
        if interval <= 0:
            interval = self.DEFAULTINTERVAL

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
                    x_axis.append(round(x_axis[-1] + interval, self.DEFAULTINTERVALROUND))

                    if defineLim:
                        maxVals, minVals = self.checkMinMax(matrixVariables, maxVals, minVals)

                    result = ebeam.getXYZ(matrixVariables)
                    twiss = result[3]
                    plot6dValues.update({x_axis[-1]: result})

                    # Aggregate the beam properties results together in a single pandas frame
                    for j, axis in enumerate(twiss.index):
                        twiss_axis = twiss.loc[axis]
                        for label, value in twiss_axis.items():
                            twiss_aggregated_df.at[axis, label].append(value)

                    # Update the progress bar
                    pbar.update(1)
                    intTrack -= interval

                if intTrack > 0:
                    # Use remainder length once when it's smaller than interval
                    matrixVariables = np.array(beamSegments[i].useMatrice(matrixVariables, length=intTrack))
                    x_axis.append(round(x_axis[-1] + intTrack, self.DEFAULTINTERVALROUND))

                    if defineLim:
                        maxVals, minVals = self.checkMinMax(matrixVariables, maxVals, minVals)

                    result = ebeam.getXYZ(matrixVariables)
                    twiss = result[3]
                    plot6dValues.update({x_axis[-1]: result})
                    # Aggregate the beam properties results together in a single pandas frame
                    for j, axis in enumerate(twiss.index):
                        twiss_axis = twiss.loc[axis]
                        for label, value in twiss_axis.items():
                            twiss_aggregated_df.at[axis, label].append(value)

                    # Update the progress bar for the remaining part
                    pbar.update(1)

        
        if saveData:
            #  Optionally save standard deviation and mean data
            name = "simulator-data-" + datetime.datetime.now().strftime('%Y-%m-%d') + "_" + datetime.datetime.now().strftime('%H_%M_%S') +".csv"
            self._csvWriteData(name, twiss_aggregated_df, x_axis)
            
        #  Testing purposes
        self.matrixVariables = matrixVariables
        self.sixdValues = plot6dValues
        if matchScaling and defineLim:
            self.setEqualAxisScaling(maxVals, minVals)

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
            # for item in twiss_aggregated_df:print(twiss_aggregated_df[item])
            colors = ['dodgerblue', 'crimson','yellow','green']
            for i in range(0,2):
                axis = twiss_aggregated_df.index[i]
                emittance = (10 ** -6) * np.array(twiss_aggregated_df.at[axis, twiss_aggregated_df.keys()[0]])
                beta = np.array(twiss_aggregated_df.at[axis, twiss_aggregated_df.keys()[2]])
                envelope = (10 ** 3) * np.sqrt(emittance * beta)
                ax5.plot(x_axis, envelope,
                         color=colors[i], linestyle='-',
                         label=r'$E_' + axis + '$ (mm)')
                


            ax5.set_ylabel('Dispersion $D$ (m)')
            ax5.set_xticks(x_axis)
            plt.xlim(0,x_axis[-1])

            #  Auto space x tick labels for readibility
            if spacing:
                totalLen = x_axis[-1]
                lastTick = x_axis[0]
                xTickLab = [lastTick]
                for tick in x_axis[1:]: 
                    if (tick - lastTick)/totalLen > self.DEFAULTSPACINGPERCENTAGE:
                        xTickLab.append(tick)
                        lastTick = tick
                    else:
                        xTickLab.append("")
                ax5.set_xticklabels(xTickLab,rotation=45,ha='right')
            else:
                ax5.set_xticklabels(x_axis,rotation=45,ha='right')

            plt.tick_params(labelsize = 9)
            plt.xlabel("Distance from start of beam (m)")
            plt.ylabel("Envelope $E$ (mm)")
            ax5.legend(loc='upper left')
            line = None
            lineList = []
            ax6 = ax5.twinx()
            for i in range(0, 2):
                axis = twiss_aggregated_df.index[i]
                dispersion = np.array(twiss_aggregated_df.at[axis, twiss_aggregated_df.keys()[4]])
                line, = ax6.plot(x_axis, dispersion,
                             color=colors[i], linestyle='--',
                             label=r'$D_' + axis + '$ (mm)')
                lineList.append(line)
                
             #for i in range(0, 4):
            #    dispersion = np.array(twiss_aggregated_df.at[twiss_aggregated_df.index[2], twiss_aggregated_df.keys()[i]])
            #    ax5.plot(x_axis, dispersion,
            #                 color=colors[i], linestyle='-',
            #                 label=twiss_aggregated_df.keys()[i])
            ax6.set_ylabel('Dispersion $D$ (mm)')

            ax6.legend(loc='upper right')

            #  Create visual representation of beamline segments
            ymin, ymax = ax5.get_ylim()
            ax5.set_ylim(ymin-(ymax*0.05), ymax)
            ymin, ymax = ax5.get_ylim()
            blockstart = 0
            for seg in beamSegments:
                rectangle = patches.Rectangle((blockstart, ymin), seg.length, ymax*0.05, linewidth=1, edgecolor=seg.color, facecolor= seg.color)
                ax5.add_patch(rectangle)
                blockstart += seg.length

            #  Important to leave tight_layout before scrollbar creation
            plt.suptitle("Beamline Simulation")
            plt.tight_layout()
            
            #   Scroll bar creation and function
            dimensions = ax5.get_position().bounds
            scrollax = plt.axes([dimensions[0],0.01,dimensions[2],0.01], facecolor = 'lightgoldenrodyellow')
            scrollbar = Slider(scrollax, 'z: 0', 0, x_axis[-1], valinit = 0, valstep=np.array(x_axis))
            scrollbar.valtext.set_visible(False)
            def update_scroll(val):
                matrix = plot6dValues.get(scrollbar.val)
                ax1.clear()
                ax2.clear()
                ax3.clear()
                ax4.clear()
                scrollbar.label.set_text("z: " + str(val))
                ebeam.plotXYZ(matrix[2], matrix[0], matrix[1], matrix[3], ax1,ax2,ax3,ax4, maxVals, minVals, defineLim, shape)
                fig.canvas.draw_idle()
            scrollbar.on_changed(update_scroll)

            lineTwissData = []
            twissDataNames = []
            for key in twiss_aggregated_df.keys():
                    lineTwissData.append([twiss_aggregated_df.at['x', key],twiss_aggregated_df.at['y', key]])
                    twissDataNames.append(key)

            class CircularList:
                index = 0

                def drawNewLines(self, ind):
                    data = lineTwissData[ind%len(lineTwissData)]
                    label = twissDataNames[ind%len(twissDataNames)]
                    for i, axis in enumerate(['x', 'y']): 
                        line = lineList[i]
                        line.set_ydata(data[i])

                        # NOTE: the legend formatting below may be changed in the future if the twiss data frame  name is changed.
                        line.set_label(label.split(' ')[0] + '$_' + axis + '$')

                    ax6.relim()
                    ax6.autoscale_view()
                    ax6.set_ylabel(label)  #  WORK AND IMPROVE THIS
                    ax6.legend(loc='upper right')
                    plt.draw()

                def nextL(self, event):
                    self.index = self.index + 1
                    self.drawNewLines(self.index)

                def prevL(self, event):
                    self.index = self.index - 1
                    self.drawNewLines(self.index)

            axprev = fig.add_axes([dimensions[0]+dimensions[2]+0.02, dimensions[1]-0.04, 0.03, 0.03])
            axnext = fig.add_axes([dimensions[0]+dimensions[2]+0.05, dimensions[1]-0.04, 0.03, 0.03])

            # axprev = fig.add_axes([])
            bnext = Button(axnext, 'Next')
            bprev = Button(axprev, 'Prev')
            circList = CircularList()
            bnext.on_clicked(circList.nextL)
            bprev.on_clicked(circList.prevL)
            
            #  Ploting
            plt.show()
        return twiss_aggregated_df