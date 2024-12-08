#   Authors: Christian Komo, Niels Bidault
from sympy import symbols, Matrix
import sympy as sp
import numpy as np

#NOTE: getSymbolicMatrice() must use all sympy methods and functions, NOT numpy
class beamline:
    def __init__(self):
        self.C = 299792458.0  # Speed of light in vacuum (m/s)
        self.q_e = 1.602176634e-19  # Elementary charge (C)
        self.m_e = 9.1093837139e-31  # Electron Mass (kg)
        self.m_p = 1.67262192595e-27  # Proton Mass (kg)
        self.m_amu = 1.66053906892E-27  # Atomic mass unit (kg)
        
        self.k_MeV = 1e-6 / self.q_e  # Conversion factor (MeV / J)

        #  [Mass (kg), charge (C), rest energy (MeV)]
        self.PARTICLES = {"electron": [self.m_e, self.q_e, (self.m_e * self.C ** 2) * self.k_MeV],
                          "proton": [self.m_p, self.q_e, (self.m_p * self.C ** 2) * self.k_MeV]}

    def changeBeamType(self, beamSegments, particleType, kineticE):
        newBeamline = beamSegments
        try:
            particleData = self.PARTICLES[particleType]
            for seg in newBeamline:
                seg.setMQE(particleData[0],particleData[1],particleData[2])
                seg.setE(kineticE)
            return newBeamline
        except KeyError:
            #  Try look for isotope particle format, format = "(isotope number),(ion charge)"
            #  ex. C12+5 (carbon 12, 5+ charge) = "12,5"
            try:
                isotopeData = particleType.split(",")
                A = int(isotopeData[0])
                Z = int(isotopeData[1])
                m_i = A * self.m_amu
                q_i = Z * self.q_e
                meV = (m_i * self.C ** 2) * self.k_MeV

                for seg in newBeamline:
                    seg.setMQE(m_i, q_i, meV)
                    seg.setE(kineticE)
                return newBeamline
            except:
                raise TypeError("Invalid particle type/isotope")


class lattice:
    #  by default every beam type is an electron beam type
    def __init__(self, length, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45):
        '''
        parent class for beamline segment object

        Parameters
        ----------
        length: float
            sets length of beamline element
        E: float, optional
            Kinetic energy value (MeV/c^2)
        '''
        self.E = E  # Kinetic energy (MeV/c^2)
        ## this should be passed from ebeam
        self.E0 = E0 # Electron rest energy (MeV/c^2)
        self.Q = Q  # (C)
        self.M = M  # (kg)
        self.C = 299792458  # Celerity (m/s)
        self.f = 2856 * (10 ** 6)  # RF frequency (Hz)
        self.gamma = (1 + (self.E / self.E0))
        self.beta = np.sqrt(1 - (1 / (self.gamma ** 2)))
        self.unitsF = 10 ** 6 # Units factor used for conversions from (keV) to (ns)
        ##
        self.color = 'none'  #Color of beamline element when graphed

        if not length <= 0:
            self.length = length
        else:
            raise ValueError("Invalid Parameter: Please enter a positive length parameter")
        
    def setE(self, E):
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

    def setMQE(self, mass, charge, restE):
        self.M = mass
        self.Q = charge
        self.E0 = restE
        self.gamma = (1 + (self.E/self.E0))
        self.beta = np.sqrt(1-(1/(self.gamma**2)))

    def getSymbolicMatrice(self):
        raise NameError("getSymbolicMatrice not defined in child class")

    def useMatrice(self, values, matrice = None):
        ''''
        Simulates the movement of given particles through its child segment

        Parameters
        ----------
        values: np.array(list[float][float])
            2D numPy list of particle elements from which to simulate data
        matrice: np.array(list[float][float])
            6x6 matrice of child segment to simulate particle data with
        '''
        if (matrice is None): 
            raise NameError("useMatrice not defined in child class")
        newMatrix = []
        for array in values:
            tempArray = np.matmul(matrice, array)
            newMatrix.append(tempArray.tolist())
        return newMatrix #  return 2d list
    
class driftLattice(lattice):
    def __init__(self, length: float, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45):
        '''
        drift lattice segment

        length: float
            length of drift segment
        '''
        super().__init__(length, E0, Q, M, E)
        self.color = "white"
        
    def getSymbolicMatrice(self, length = None):
        if length is None: l = self.length
        else: l = symbols(length, real = True)
        M56 = (l * self.f / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        mat = Matrix([[1, l, 0, 0, 0, 0],
                      [0, 1, 0, 0, 0, 0],
                      [0, 0, 1, l, 0, 0],
                      [0, 0, 0, 1, 0, 0],
                      [0, 0, 0, 0, 1, M56],
                      [0, 0, 0, 0, 0, 1]])
        return mat
    
    #Matrix multiplecation, values is a 2 dimensional numPy array, each array is 6 elements long
    #values = np.array([[x, x', y, y', z, z'],...])
    #Note: 1x6 array is multiplied correctly with 6x6 array
    def useMatrice(self, values, length = -1):
        if length <= 0:
            length = self.length
        M56 = (length * self.f / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        return super().useMatrice(values,(np.array([[1, length, 0, 0, 0, 0],
                                                    [0, 1, 0, 0, 0, 0],
                                                    [0, 0, 1, length, 0, 0],
                                                    [0, 0, 0, 1, 0, 0],
                                                    [0, 0, 0, 0, 1, M56],
                                                    [0, 0, 0, 0, 0, 1]])))

    def testSymbol(self, length = None, numeric = False):
        mat = Matrix([[1, l, 0, 0, 0, 0],
                      [0, 1, 0, 0, 0, 0],
                      [0, 0, 1, l, 0, 0],
                      [0, 0, 0, 1, 0, 0],
                      [0, 0, 0, 0, 1, M56],
                      [0, 0, 0, 0, 0, 1]])

        if numeric:
            return np.array([[1, length, 0, 0, 0, 0],
                                                    [0, 1, 0, 0, 0, 0],
                                                    [0, 0, 1, length, 0, 0],
                                                    [0, 0, 0, 1, 0, 0],
                                                    [0, 0, 0, 0, 1, M56],
                                                    [0, 0, 0, 0, 0, 1]])

        if length is None: l = self.length
        else: l = symbols(length, real = True)
        M56 = (l * self.f / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        mat = Matrix([[1, l, 0, 0, 0, 0],
                      [0, 1, 0, 0, 0, 0],
                      [0, 0, 1, l, 0, 0],
                      [0, 0, 0, 1, 0, 0],
                      [0, 0, 0, 0, 1, M56],
                      [0, 0, 0, 0, 0, 1]])
        return mat

    def __str__(self):
        return f"Drift beamline segment {self.length} m long"


class qpfLattice(lattice):
    def __init__(self, current: float, length: float = 0.0889, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45):
        super().__init__(length, E0, Q, M, E)
        self.current = current # Amps
        self.color = "cornflowerblue"
        self.G = 2.694  # Quadruple focusing strength (T/A/m)

    def getSymbolicMatrice(self, length = None, current = None):
        if current is None: I = self.current
        else: I = symbols(current, real = True)
        if length is None: l = self.length
        else: l = symbols(length, real = True)
        

        self.k = sp.Abs((self.Q*self.G*I)/(self.M*self.C*self.beta*self.gamma))
        self.theta = sp.sqrt(self.k)*l

        M11 = sp.cos(self.theta)
        M21 = (-(sp.sqrt(self.k)))*sp.sin(self.theta)
        M22 = sp.cos(self.theta)
        M33 = sp.cosh(self.theta)
        M43 = sp.sqrt(self.k)*sp.sinh(self.theta)
        M44 = sp.cosh(self.theta)
        M56 = (l * self.f / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        if I == 0:
            M12 = l
            M34 = l
        else:
            M34 = sp.sinh(self.theta)*(1/sp.sqrt(self.k))
            M12 = sp.sin(self.theta)*(1/sp.sqrt(self.k))

        mat =  Matrix([[M11, M12, 0, 0, 0, 0],
                        [M21, M22, 0, 0, 0, 0],
                        [0, 0, M33, M34, 0, 0],
                        [0, 0, M43, M44, 0, 0],
                        [0, 0, 0, 0, 1, M56],
                        [0, 0, 0, 0, 0, 1]])
        
        return mat


    '''
    performs a transformation to a 2d np array made of 1x6 variable matrices

    values: np.array([list[int],...])
    '''
    def useMatrice(self, values, length = -1, current = -1):
        if length <= 0:
            length = self.length
        if current < 0:
            current = self.current

        #   Necessary because code had problems working with numpy arrays
        if isinstance(current, np.ndarray):
            current = current[0]

        self.k = np.abs((self.Q*self.G*current)/(self.M*self.C*self.beta*self.gamma))
        self.theta = np.sqrt(self.k)*length

        M11 = np.cos(self.theta)
        M21 = (-(np.sqrt(self.k)))*np.sin(self.theta)
        M22 = np.cos(self.theta)
        M33 = np.cosh(self.theta)
        M43 = np.sqrt(self.k)*np.sinh(self.theta)
        M44 = np.cosh(self.theta)
        M56 = (length * self.f / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        if current == 0:
            M12 = length
            M34 = length
        else:
            M34 = np.sinh(self.theta) * (1 / np.sqrt(self.k))
            M12 = np.sin(self.theta) * (1 / np.sqrt(self.k))

        return super().useMatrice(values, np.array([[M11, M12, 0, 0, 0, 0],
                                                    [M21, M22, 0, 0, 0, 0],
                                                    [0, 0, M33, M34, 0, 0],
                                                    [0, 0, M43, M44, 0, 0],
                                                    [0, 0, 0, 0, 1, M56],
                                                    [0, 0, 0, 0, 0, 1]]))
    
    def __str__(self):
        return f"QPF beamline segment {self.length} mm long"


class qpdLattice(lattice):
    def __init__(self, current: float, length: float = 0.0889, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45):
        super().__init__(length, E0, Q, M, E)
        self.current = current # Amps
        self.G = 2.694  # Quadruple focusing strength (T/A/m)
        self.color = "lightcoral"
        
    def getSymbolicMatrice(self, length = None, current = None):
        if current is None: I = self.current
        else: I = symbols(current, real = True)
        if length is None: l = self.length
        else: l = symbols(length, real = True)
        
        self.k = sp.Abs((self.Q*self.G*I)/(self.M*self.C*self.beta*self.gamma))
        self.theta = sp.sqrt(self.k)*l

        M11 = sp.cosh(self.theta)
        M21 = sp.sqrt(self.k)*sp.sinh(self.theta)
        M22 = sp.cosh(self.theta)
        M33 = sp.cos(self.theta)
        M43 = (-(sp.sqrt(self.k)))*sp.sin(self.theta)
        M44 = sp.cos(self.theta)
        M56 = l * self.f / (self.C * self.beta * self.gamma * (self.gamma + 1))

        if I == 0:
            M12 = l
            M34 = l
        else:
            M34 = sp.sin(self.theta)*(1/sp.sqrt(self.k))
            M12 = sp.sinh(self.theta)*(1/sp.sqrt(self.k))

        mat = Matrix([[M11, M12, 0, 0, 0, 0],
                        [M21, M22, 0, 0, 0, 0],
                        [0, 0, M33, M34, 0, 0],
                        [0, 0, M43, M44, 0, 0],
                        [0, 0, 0, 0, 1, M56],
                        [0, 0, 0, 0, 0, 1]])
        
        return mat

    def useMatrice(self, values, length = -1, current = -1):
        if length <= 0:
            length = self.length
        if current < 0:
            current = self.current

        #   Necessary because code had problems working with numpy arrays
        if isinstance(current, np.ndarray):
            current = current[0]

        self.k = np.abs((self.Q*self.G*current)/(self.M*self.C*self.beta*self.gamma))
        self.theta = np.sqrt(self.k)*length

        M11 = np.cosh(self.theta)
        M21 = np.sqrt(self.k)*np.sinh(self.theta)
        M22 = np.cosh(self.theta)
        M33 = np.cos(self.theta)
        M43 = (-(np.sqrt(self.k)))*np.sin(self.theta)
        M44 = np.cos(self.theta)
        M56 = length * self.f / (self.C * self.beta * self.gamma * (self.gamma + 1))
        if current == 0:
            M12 = length
            M34 = length
        else:
            M34 = np.sin(self.theta) * (1 / np.sqrt(self.k))
            M12 = np.sinh(self.theta) * (1 / np.sqrt(self.k))

        return super().useMatrice(values, (np.array([[M11, M12, 0, 0, 0, 0],
                                                    [M21, M22, 0, 0, 0, 0],
                                                    [0, 0, M33, M34, 0, 0],
                                                    [0, 0, M43, M44, 0, 0],
                                                    [0, 0, 0, 0, 1, M56],
                                                    [0, 0, 0, 0, 0, 1]])))

    def __str__(self):
        return f"QPD beamline segment {self.length} m long"

class dipole(lattice):
    def __init__(self, length: float = 0.0889, angle: float = 1.5, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45):
        super().__init__(length, E0, Q, M, E)
        self.color = "forestgreen"
        self.angle = angle  # degrees
        self.length_total = length

    def getSymbolicMatrice(self, length = None, angle = None):
        if length is None: L = self.length
        else: L = symbols(length, real = True)
        if angle is None: a = self.angle
        else: a = symbols(angle, real = True)

        by = (self.M*self.C*self.beta*self.gamma / self.Q) * (a * sp.pi / 180 / self.length_total)
        rho = self.M*self.C*self.beta*self.gamma / (self.Q * by)
        theta = L / rho
        C = sp.cos(theta)
        S = sp.sin(theta)

        M16 = rho * (1 - C) * (self.gamma / (self.gamma + 1))
        M26 = S * (self.gamma / (self.gamma + 1))
        M51 = self.f * (-S / (self.beta * self.C))
        M52 = self.f * (-rho * (1 - C) / (self.beta * self.C))
        M56 = self.f * (-rho * (L / rho - S) / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        mat = Matrix([[C, rho * S, 0, 0, 0, M16],
                      [-S / rho, C, 0, 0, 0, M26],
                      [0, 0, 1, L, 0, 0],
                      [0, 0, 0, 1, 0, 0],
                      [M51, M52, 0, 0, 1, M56],
                      [0, 0, 0, 0, 0, 1]])
        
        return mat

    
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
        By = (self.M*self.C*self.beta*self.gamma / self.Q) * (angle * np.pi / 180 / self.length_total)

        rho = self.M*self.C*self.beta*self.gamma / (self.Q * By)

        theta = length / rho
        C = np.cos(theta)
        S = np.sin(theta)
        L = length

        M16 = rho * (1 - C) * (self.gamma / (self.gamma + 1))
        M26 = S * (self.gamma / (self.gamma + 1))
        M51 = self.f * (-S / (self.beta * self.C))
        M52 = self.f * (-rho * (1 - C) / (self.beta * self.C))
        M56 = self.f * (-rho * (L / rho - S) / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        M = np.array([[C, rho * S, 0, 0, 0, M16],
                      [-S / rho, C, 0, 0, 0, M26],
                      [0, 0, 1, length, 0, 0],
                      [0, 0, 0, 1, 0, 0],
                      [M51, M52, 0, 0, 1, M56],
                      [0, 0, 0, 0, 0, 1]])

        return super().useMatrice(values, M)

    def __str__(self):
        return f"Horizontal dipole magnet segment {self.length} m long (curvature)"

class dipole_wedge(lattice):
    def __init__(self, length, angle: float = 1, dipole_length: float = 0.0889, dipole_angle: float = 1.5, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45):
        super().__init__(length, E0, Q, M, E)
        self.color = "lightgreen"
        self.angle = angle
        self.dipole_length = dipole_length
        self.dipole_angle = dipole_angle
        self.gap = 0.01  # hard-edge fringe field gap (m)


    def getSymbolicMatrice(self, length = None, angle = None):
        if length is None: L = self.length
        else: L = symbols(length, real = True)
        if angle is None: a = self.angle
        else: a = symbols(angle, real = True)
        dipole_angle = self.dipole_angle
        dipole_length = self.dipole_length

        # Hard edge model for the wedge magnets
        By = (self.M*self.C*self.beta*self.gamma / self.Q) * (dipole_angle * sp.pi / 180 / dipole_length)
        R = self.M*self.C*self.beta*self.gamma / (self.Q * By)
        k = sp.Abs((self.Q * By / self.length) / (self.M * self.C * self.beta * self.gamma)) # Verify
        eta = (a * sp.pi / 180) * L / self.length # Verify
        E = (L * k) / ((sp.cos(eta)) ** 2)
        T = sp.tan(eta)
        M56 = self.f * (L / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        mat = Matrix([[1, 0, 0, 0, 0, 0],
                      [-T / R, 1, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, -(T + E / 3) / R, 1, 0, 0],
                      [0, 0, 0, 0, 1, M56],
                      [0, 0, 0, 0, 0, 1]])
        
        return mat

    '''
    performs a transformation to a 2d np array made of 1x6 variable matrices

    values: np.array([list[int],...])
    '''
    def useMatrice(self, values, length = -1, angle = -1):
        if length <= 0:
            length = self.length
        if angle < 0:
            angle = self.angle
        dipole_angle = self.dipole_angle
        dipole_length = self.dipole_length

        #   TEMPORARY PURPOSES
        if isinstance(angle, np.ndarray):
            angle = angle[0]

        # Hard edge model for the wedge magnets
        By = (self.M*self.C*self.beta*self.gamma / self.Q) * (dipole_angle * np.pi / 180 / dipole_length)
        R = self.M*self.C*self.beta*self.gamma / (self.Q * By)
        k = np.abs((self.Q * By / self.length) / (self.M * self.C * self.beta * self.gamma)) # Verify
        eta = (angle * np.pi / 180) * length / self.length # Verify
        E = (length * k) / ((np.cos(eta)) ** 2)
        T = np.tan(eta)
        M56 = self.f * (length / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        M = np.array([[1, 0, 0, 0, 0, 0],
                      [-T / R, 1, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, -(T + E / 3) / R, 1, 0, 0],
                      [0, 0, 0, 0, 1, M56],
                      [0, 0, 0, 0, 0, 1]])

        return super().useMatrice(values, M)
    def __str__(self):
        return f"Horizontal wedge dipole magnet segment {self.length} m long (curvature)"
    
