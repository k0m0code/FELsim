#   Authors: Christian Komo, Niels Bidault

import numpy as np

#should I use the position that each beamline object takes up instead of its total length?
#Give each beamline object its optional plot6d paramter at certain lengths into the beamline?

class lattice:
    def __init__(self, length: float = 0, E = 35):
        '''
        parent class for beamline segment object

        Parameters
        ----------
        length: float
            sets length of beamline element
        E: float, optional
            Kinetic energy value (???)
        '''
        self.E = E  # Kinetic energy (MeV/c^2)
        self.E0 = 0.51099
        self.QE = 1.60217663e-19  #C
        self.ME = 9.1093837e-31  #kg
        self.C = 299792458  #m/s
        self.G = 2.694  #Quadruple focusing strength (T/A/m)
        self.gamma = (1 + (self.E/self.E0))
        self.beta = np.sqrt(1 - (1 / (self.gamma ** 2)))
        self.color = 'none'  #Color of beamline element when graphed
        if not length <= 0:
            self.length = length
        else:
            raise ValueError("Invalid Parameter: Please enter a positive length parameter")
        
    def changeE(self, E):
        '''
        Set value of E and its subsequent dependent formulas

        Parameters
        ----------
        E: float
            New value of E
        '''
        self.E = E
        self.gamma = (1 + (self.E/self.E0))
        self.beta = np.sqrt(1-(1/(self.gamma**2)))

    def useMatrice(self, values, matrice):
        ''''
        Simulates the movement of given particles through its child segment

        Parameters
        ----------
        values: np.array(list[float][float])
            2D numPy list of particle elements from which to simulate data
        matrice: np.array(list[float][float])
            6x6 matrice of child segment to simulate particle data with
        '''
        newMatrix = []
        for array in values:
            tempArray = np.matmul(matrice, array)
            newMatrix.append(tempArray.tolist())
        return newMatrix #  return 2d list

class driftLattice(lattice):
    def __init__(self, length: float = -1):
        '''
        drift lattice segment

        length: float
            length of drift segment
        '''
        super().__init__(length = length)
        self.color = "white"
        
    #Matrix multiplecation, values is a 2 dimensional numPy array, each array is 6 elements long
    #values = np.array([[x, x', y, y', z, z'],...])
    #Note: 1x6 array is multiplied correctly with 6x6 array
    def useMatrice(self, values, length = -1):
        if length <= 0:
            length = self.length
        return super().useMatrice(values,(np.array([[1, length, 0, 0, 0, 0],
                                 [0, 1, 0, 0, 0, 0],
                                 [0, 0, 1, length, 0, 0],
                                 [0, 0, 0, 1, 0, 0],
                                 [0, 0, 0, 0, 1, (length/((self.gamma)**2))],
                                 [0, 0, 0, 0, 0, 1]])))
    
    def __str__(self):
        return f"Drift beamline segment {self.length} m long"


class qpfLattice(lattice):
    def __init__(self, current: float, length: float = 0.0889):
        super().__init__(length = length)
        self.current = current # Amps
        self.color = "cornflowerblue"
    '''
    performs a transformation to a 2d np array made of 1x6 variable matrices

    values: np.array([list[int],...])
    '''
    def useMatrice(self, values, length = -1, current = -1):
        if length <= 0:
            length = self.length
        if current < 0:
            current = self.current

        #   TEMPORARY PURPOSES
        if isinstance(current, np.ndarray):
            current = current[0]

        self.k = np.abs((self.QE*self.G*current)/(self.ME*self.C*self.beta*self.gamma))
        self.theta = np.sqrt(self.k)*length

        field1 = np.cos(self.theta)
        field2 = np.sin(self.theta)*(1/np.sqrt(self.k))
        field3 = (-(np.sqrt(self.k)))*np.sin(self.theta)
        field4 = np.cos(self.theta)
        field5 = np.cosh(self.theta)
        field6 = np.sinh(self.theta)*(1/np.sqrt(self.k))
        field7 = np.sqrt(self.k)*np.sinh(self.theta)
        field8 = np.cosh(self.theta)
        field9 = length/(self.gamma**2)

        return super().useMatrice(values, np.array([[field1, field2, 0, 0, 0, 0],
                                                    [field3, field4, 0, 0, 0, 0],
                                                    [0, 0, field5, field6, 0, 0],
                                                    [0, 0, field7, field8, 0, 0],
                                                    [0, 0, 0, 0, 1, field9],
                                                    [0, 0, 0, 0, 0, 1]]))
    
    def __str__(self):
        return f"QPF beamline segment {self.length} mm long"


class qpdLattice(lattice):
    def __init__(self, length: float = 0.0889, current: float = 0):
        super().__init__(length)
        self.current = current # Amps
        self.color = "lightcoral"

    def useMatrice(self, values, length = -1, current = -1):
        if length <= 0:
            length = self.length
        if current < 0:
            current = self.current

        #   TEMPORARY PURPOSES
        if isinstance(current, np.ndarray):
            current = current[0]

        self.k = np.abs((self.QE*self.G*current)/(self.ME*self.C*self.beta*self.gamma))
        self.theta = np.sqrt(self.k)*length

        field1 = np.cos(self.theta)
        field2 = np.sin(self.theta)*(1/np.sqrt(self.k))
        field3 = (-(np.sqrt(self.k)))*np.sin(self.theta)
        field4 = np.cos(self.theta)
        field5 = np.cosh(self.theta)
        field6 = np.sinh(self.theta)*(1/np.sqrt(self.k))
        field7 = np.sqrt(self.k)*np.sinh(self.theta)
        field8 = np.cosh(self.theta)
        field9 = length/(self.gamma**2)

        return super().useMatrice(values, (np.array([[field5, field6, 0, 0, 0, 0],
                                                    [field7, field8, 0, 0, 0, 0],
                                                    [0, 0, field1, field2, 0, 0],
                                                    [0, 0, field3, field4, 0, 0],
                                                    [0, 0, 0, 0, 1, field9],
                                                    [0, 0, 0, 0, 0, 1]])))

    def __str__(self):
        return f"QPD beamline segment {self.length} m long"

class dipole(lattice):
    def __init__(self, length: float = 0.0889, angle: float = 1.5):
        super().__init__(length=length)
        self.color = "forestgreen"
        self.angle = angle # degrees
        self.By = (self.ME*self.C*self.beta*self.gamma / self.QE) * (self.angle * np.pi / 180 / self.length)
        self.rho = self.ME*self.C*self.beta*self.gamma / (self.QE * self.By)

    '''
    performs a transformation to a 2d np array made of 1x6 variable matrices

    values: np.array([list[int],...])
    '''
    def useMatrice(self, values, length = -1, angle = -1):
        if length <= 0:
            length = self.length
        if angle < 0:
            angle = self.angle

        #   TEMPORARY PURPOSES
        if isinstance(angle, np.ndarray):
            angle = angle[0]

        # Rectangular dipole
        R = self.rho
        theta = length / self.rho
        C = np.cos(theta)
        S = np.sin(theta)
        L = length
        G = self.gamma

        M = np.array([[C, R * S, 0, 0, 0, R * (1 - C)],
                      [-S / R, C, 0, 0, 0, S],
                      [0, 0, 1, L, 0, 0],
                      [0, 0, 0, 1, 0, 0],
                      [S, R * (1 - C), 0, 0, 1, R * (L / R - S) + L / (G ** 2)],
                      [0, 0, 0, 0, 0, 1]])

        return super().useMatrice(values, M)

    def __str__(self):
        return f"Horizontal dipole magnet segment {self.length} m long (curvature)"

class dipole_wedge(lattice):
    def __init__(self, length, angle: float = 1, dipole_length: float = 0.0889, dipole_angle: float = 1.5):
        super().__init__(length=length)
        self.color = "lightgreen"
        self.angle = angle
        self.By = (self.ME*self.C*self.beta*self.gamma / self.QE) * (dipole_angle * np.pi / 180 / dipole_length)
        self.rho = self.ME*self.C*self.beta*self.gamma / (self.QE * self.By)
        self.gap = 0.01  # hard-edge fringe field gap (m)

    '''
    performs a transformation to a 2d np array made of 1x6 variable matrices

    values: np.array([list[int],...])
    '''
    def useMatrice(self, values, length = -1, angle = -1):
        if length <= 0:
            length = self.length
        if angle < 0:
            angle = self.angle

        #   TEMPORARY PURPOSES
        if isinstance(angle, np.ndarray):
            angle = angle[0]

        # Hard edge model for the wedge magnets
        R = self.rho
        k = np.abs((self.QE * self.By / self.length) / (self.ME * self.C * self.beta * self.gamma)) # Verify
        eta = (self.angle * np.pi / 180) * length / self.length # Verify
        E = (length * k) / ((np.cos(eta)) ** 2)
        T = np.tan(eta)
        L = length
        G = self.gamma

        M = np.array([[1, 0, 0, 0, 0, 0],
                      [-T / R, 1, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, -(T + E / 3) / R, 1, 0, 0],
                      [0, 0, 0, 0, 1, L / G ** 2],
                      [0, 0, 0, 0, 0, 1]])


        return super().useMatrice(values, M)
    def __str__(self):
        return f"Horizontal wedge dipole magnet segment {self.length} m long (curvature)"
