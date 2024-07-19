#   Authors: Niels Bidault, Christian Komo
import numpy as np
import matplotlib.pyplot as plt

#in plotDriftTransform, add legend and gausian distribution for x and y points
#Replace list variables that are unchaning with tuples, more efficient for calculations
#GetDriftMatrice should handle looping through all the diff values in the variable list of each point

class beam:
    def __init__(self, length: float):
        self.E = 35  # Kinetic energy (MeV/c^2)
        self.E0 = 0.510999
        self.length = length

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

    #Matrix multiplecation, values is a 1 dimensional numPy array or normal list 6 units long
    #values = [x, x', y, y', z, z']
    def getDriftMatrice(self, values: list[int], length = -1):
        if length != -1:
            length = self.length
        gamma = (1 + (self.E/self.E0))

        row1 = values[0] + values[1]*self.length
        row2 = values[1]
        row3 = values[2] + values[3]*self.length
        row4 = values[3]
        row5 = values[4] + (values[5]*(self.length/(gamma**2)))
        row6 = values[5]
        newMatrix = [row1,row2,row3,row4,row5,row6]
        return newMatrix
    
    # Can length variable be negative?
    def DriftTransform(self, x_pos: list[int], y_pos: list[int], phase_x: list[int], phase_y: list[int], length = -1, plot = True):
        if length == -1:
            length = self.length
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
    def plotBeamPositionTransform(self, matrixVariables, beamSegments, interval):
        x_pos = matrixVariables[:, 0]
        y_pos = matrixVariables[:, 2]
        phase_x = matrixVariables[:, 1]
        phase_y = matrixVariables[:, 3]
        xUpdated = [np.std(x_pos)]
        yUpdated = [np.std(y_pos)]
        x_axis = [0]

        for i in range(len(beamSegments[0])):
            if beamSegments[0][i] == "drift":
                intTrack = beamSegments[1][i]
                while intTrack > interval:
                    xtrans, ytrans = self.DriftTransform(x_pos,y_pos,phase_x,phase_y, length = interval,plot = False)
                    xUpdated.append(np.std(xtrans))
                    yUpdated.append(np.std(ytrans))
                    x_pos = xtrans
                    y_pos = ytrans
                    intTrack = intTrack - interval
                    x_axis.append(x_axis[len(x_axis)-1]+interval)
                xtrans, ytrans = self.DriftTransform(x_pos,y_pos,phase_x,phase_y, length = intTrack,plot = False)
                xUpdated.append(np.std(xtrans))
                yUpdated.append(np.std(ytrans))
                x_pos = xtrans
                y_pos = ytrans
                x_axis.append(x_axis[len(x_axis)-1]+intTrack)
            #Add more matrix multiplcation down here
        
        fig, ax = plt.subplots()
        plt.plot(x_axis, xUpdated)
        plt.plot(x_axis, yUpdated)
        ax.set_xticks(x_axis)
        plt.xlim(0,x_axis[len(x_axis)-1] + x_axis[len(x_axis)-1]*0.10)
        ax.set_xticklabels(x_axis,rotation=45,ha='right')
        plt.tick_params(labelsize = 9)
        plt.xlabel("Distance from start of beam (mm??)")    # Change the units of xlabel?
        plt.ylabel("Standard deviation (mm??)") # Change the units of ylabel?
        plt.show()        

