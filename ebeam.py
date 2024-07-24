#   Authors: Niels Bidault, Christian Komo
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import csv
import datetime

#in plotDriftTransform, add legend and gausian distribution for x and y points
#Replace list variables that are unchanging with tuples, more efficient for calculations
#For functions like drifttransform, paramter should only be a single 2d array with all 6 initial variables
#GetDriftMatrice should handle looping through all the diff values in the variable list of each point
#Add legend for graphs like plotBeamPositionTransform
#Replace manual multiplecation with matrice multiplcation
#Use getDriftMatrice instead of drifttransform in plotbeampoisitiontransform
#Same total pipe length but different interval should end with the same standard deviation
#EACH BEAMLINE OBJECT SHOULD REPRESENT A DIFFERENT SECTION OF THE BEAM??



class beam:
    def __init__(self, driftLength: float = 0, qpfLength: float = 88.9, current: float = 0):
        self.E = 35  # Kinetic energy (MeV/c^2)
        self.E0 = 0.51099
        self.current = current
        self.driftLength = driftLength
        self.qpfLength = qpfLength

    def csvWriteData(self, name, distance, std_x, std_y, mean_x, mean_y):
        with open(name, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            if csvfile.tell() == 0:
                csvwriter.writerow(['distance (mm)','x position standard deviaton (mm)', 'y position standard deviation (mm)', 'x position mean (mm)', 'y position mean (mm)'])
            csvwriter.writerow([distance, std_x, std_y, mean_x, mean_y])

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
    def plot_6d(self, position_x_values: list[int], position_y_values: list[int], phase_x_values: list[int], phase_y_values: list[int],
                energy_values: list[int], time_values: list[int]):
        
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

    #Matrix multiplecation, values is a 2 dimensional numPy array, each array is 6 elements long
    #values = np.array([[x, x', y, y', z, z'],...])
    #Note: 1x6 array is multiplied correctly with 6x6 array
    def getDriftMatrice(self, values, length = -1):
        if length == -1:
            length = self.driftLength
        gamma = (1 + (self.E/self.E0))

        driftMatrice = (np.array([[1, length, 0, 0, 0, 0],
                                 [0, 1, 0, 0, 0, 0],
                                 [0, 0, 1, length, 0, 0],
                                 [0, 0, 0, 1, 0, 0],
                                 [0, 0, 0, 0, 1, (length/(gamma**2))],
                                 [0, 0, 0, 0, 0, 1]]))

        newMatrix = []
        for array in values:
            tempArray = np.matmul(driftMatrice, array)
            newMatrix.append(tempArray.tolist())
        return newMatrix #  return 2d list
    
    '''
    performs a transformation to a 2d np array made of 1x6 variable matrices

    values: np.array([list[int],...])
    '''
    def getQPFmatrice(self, values, length = -1, current = -1):
        if length == -1:
            length = self.qpfLength
        if current == -1:
            current = self.current
        theta = np.sqrt(current)*length
        gamma = (1 + (self.E/self.E0))

        qpfMatrice = np.array([[np.cos(theta),(np.sin(theta)/np.sqrt(current)),0,0,0,0],
                               [(-(np.sqrt(current)))*(np.sin(theta)),np.cos(theta),0,0,0,0],
                               [0,0,np.cosh(theta),(np.sinh(theta))/(np.sqrt(current)),0,0],
                               [0,0,np.sqrt(current)*np.sinh(theta),np.cosh(theta),0,0],
                               [0,0,0,0,1,length/(gamma**2)],
                               [0,0,0,0,0,1]])
        
        newMatrix = []
        for array in values:
            tempArray = np.matmul(qpfMatrice, array)
            newMatrix.append(tempArray.tolist())
        return newMatrix
    
    #integrate getDriftMatrice within this function or get rid of getDriftMatrice altogether
    # Can length variable be negative?
    def driftTransformScatter(self, x_pos: list[int], y_pos: list[int], phase_x: list[int], phase_y: list[int], length = -1, plot = True):
        if length == -1:
            length = self.driftLength
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

    '''
    matrixvairables: list[float][float]
    2d numpy array containing initial condiitons

    beamSegmeents: list[str][float]
    2d numpy array or 2D list containing 1. the type of "medium" beam passes through 2. the length of each medium matching with its index
    '''
    def plotBeamPositionTransform(self, matrixVariables, beamSegments, interval = 1, saveData = False):
        xUpdated = [np.std(matrixVariables[:, 0])]
        yUpdated = [np.std(matrixVariables[:, 2])]
        xMean = [np.mean(matrixVariables[:, 0])]
        yMean = [np.mean(matrixVariables[:, 2])]
        x_axis = [0]
        xaxisMax = sum(beamSegments[1])

        for i in range(len(beamSegments[0])):
            intTrack = beamSegments[1][i]
            if beamSegments[0][i] == "drift":
                while intTrack >= interval:
                    matrixVariables = np.array(self.getDriftMatrice(matrixVariables, length = interval))
                    xUpdated.append(np.std(matrixVariables[:,0]))
                    yUpdated.append(np.std(matrixVariables[:,2]))
                    xMean.append(np.mean(matrixVariables[:,0]))
                    yMean.append(np.mean(matrixVariables[:,2]))
                    intTrack -= interval
                    x_axis.append(round(x_axis[-1]+interval, 3))
                if intTrack > 0:
                    matrixVariables = np.array(self.getDriftMatrice(matrixVariables, length = intTrack))
                    xUpdated.append(np.std(matrixVariables[:,0]))
                    yUpdated.append(np.std(matrixVariables[:,2]))
                    xMean.append(np.mean(matrixVariables[:,0]))
                    yMean.append(np.mean(matrixVariables[:,2]))
                    x_axis.append(round(x_axis[-1]+intTrack, 3))
            if beamSegments[0][i] == 'QPF':
                while intTrack >= interval:
                    matrixVariables = np.array(self.getQPFmatrice(matrixVariables, length = interval, current = 0.0001)) #test current value
                    xUpdated.append(np.std(matrixVariables[:,0]))
                    yUpdated.append(np.std(matrixVariables[:,2]))
                    xMean.append(np.mean(matrixVariables[:,0]))
                    yMean.append(np.mean(matrixVariables[:,2]))
                    intTrack -= interval
                    x_axis.append(round(x_axis[-1]+interval, 3))
                if intTrack > 0:
                    matrixVariables = np.array(self.getQPFmatrice(matrixVariables, length = intTrack, current = 0.0001)) #test current value
                    xUpdated.append(np.std(matrixVariables[:,0]))
                    yUpdated.append(np.std(matrixVariables[:,2]))
                    xMean.append(np.mean(matrixVariables[:,0]))
                    yMean.append(np.mean(matrixVariables[:,2]))
                    x_axis.append(round(x_axis[-1]+intTrack, 3))
            # Add more matrix multiplcation below

        if saveData:
            name = "simulator-data-" + datetime.datetime.now().strftime('%Y-%m-%d') + "_" + datetime.datetime.now().strftime('%H_%M_%S') +".csv"
            for i in range(len(x_axis)):
                self.csvWriteData(name, x_axis[i], xUpdated[i], yUpdated[i], xMean[i], yMean[i])

        fig, ax = plt.subplots()
        plt.plot(x_axis, xUpdated)
        plt.plot(x_axis, yUpdated)
        plt.plot(x_axis, xMean, color = 'red')
        plt.plot(x_axis, yMean, color = 'blue')
        plt.subplots_adjust(bottom=0.25)
        ax.set_xticks(x_axis)
        plt.xlim(0,x_axis[-1] + x_axis[-1]*0.10)
        ax.set_xticklabels(x_axis,rotation=45,ha='right')
        plt.tick_params(labelsize = 9)
        plt.xlabel("Distance from start of beam (mm)")
        plt.ylabel("Standard deviation (mm)")
        plt.xlim(0, xaxisMax*0.2)
        scrollax = plt.axes([0.1,0.02,0.8,0.06], facecolor = 'lightgoldenrodyellow')
        scrollbar = Slider(scrollax, 'scroll', 0, 100, valinit = 0, valstep=1)

        def update_scroll(val):
            pos = scrollbar.val
            ax.set_xlim((pos/125)*xaxisMax, (pos/125)*xaxisMax + xaxisMax*(0.2))
            fig.canvas.draw_idle()

        scrollbar.on_changed(update_scroll)
        
        plt.tight_layout()
        plt.show()
