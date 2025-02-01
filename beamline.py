#   Authors: Christian Komo, Niels Bidault
from sympy import symbols, Matrix
import sympy as sp
import numpy as np
from scipy import interpolate
import math
        
class lattice:
    #  by default every beam type is an electron beam type
    def __init__(self, length, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45, fringeType = None):
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
        self.fringeType = fringeType  # Each segment has no magnetic fringe by default
        # self.startPos = None
        # self.endPos = None

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
    
    def getSymbolicMatrice(self, **kwargs):
        raise NotImplementedError("getSymbolicMatrice not defined in child class")
    
    #unfortuately, cannot check whether the kwargs exist in the segments function or not alreay\dy
    def useMatrice(self, val, **kwargs):
        ''''
        Simulates the movement of given particles through its child segment with the 
        given numeric parameters for its segment

        Parameters
        ----------
        val: np.array(list[float][float])
            2D numPy list of particle elements from which to simulate data
        **kwargs: 
            other segment specific parameters.
        '''
        mat = self.getSymbolicMatrice(numeric = True, **kwargs)
        npMat = np.array(mat).astype(np.float64)

        newMatrix = []
        for array in val:
            tempArray = np.matmul(npMat, array)
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
    
    # note: unlike old usematrice, this func doesnt check for negative/zero parameter numbers,
    # Nor if length is actually the dtype numeric specifies
    # implement both useMatrice and symbolic matrice in this, have to delete both later
    def getSymbolicMatrice(self, numeric = False, length = None):
        l = None
        if length is None:
            l = self.length
        else:
            if numeric: l = length  # length should be number
            else: l = symbols(length, real = True)  # length should be string
        
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
    def __init__(self, current: float, length: float = 0.0889, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45, fringeType = 'decay'):
        super().__init__(length, E0, Q, M, E, fringeType)
        self.current = current # Amps
        self.color = "cornflowerblue"
        self.G = 2.694  # Quadruple focusing strength (T/A/m)

    def getSymbolicMatrice(self, numeric = False, length = None, current = None):
        l = None
        I = None

        if length is None:
            l = self.length
        else:
            if numeric: l = length  # length should be number
            else: l = symbols(length, real = True)  # length should be string

        if current is None:
            I = self.current
        else:
            if numeric: I = current  # current should be number
            else: I = symbols(current, real = True)  # current should be string

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
    
    def __str__(self):
        return f"QPF beamline segment {self.length} mm long and a current of {self.current} amps"


class qpdLattice(lattice):
    def __init__(self, current: float, length: float = 0.0889, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45, fringeType = 'decay'):
        super().__init__(length, E0, Q, M, E, fringeType)
        self.current = current # Amps
        self.G = 2.694  # Quadruple focusing strength (T/A/m)
        self.color = "lightcoral"
    
    def getSymbolicMatrice(self, numeric = False, length = None, current = None):
        l = None
        I = None

        if length is None:
            l = self.length
        else:
            if numeric: l = length  # length should be number
            else: l = symbols(length, real = True)  # length should be string

        if current is None:
            I = self.current
        else:
            if numeric: I = current  # current should be number
            else: I = symbols(current, real = True)  # current should be string

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

    def __str__(self):
        return f"QPD beamline segment {self.length} m long and a current of {self.current} amps"

class dipole(lattice):
    def __init__(self, length: float = 0.0889, angle: float = 1.5, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45, fringeType = 'decay'):
        super().__init__(length, E0, Q, M, E, fringeType)
        self.color = "forestgreen"
        self.angle = angle  # degrees
    
    def getSymbolicMatrice(self, numeric = False, length = None, angle = None):
        l = None
        a = None

        if length is None:
            l = self.length
        else:
            if numeric: l = length  # length should be number
            else: l = symbols(length, real = True)  # length should be string

        if angle is None:
            a = self.angle
        else:
            if numeric: a = angle  # angle should be number
            else: a = symbols(angle, real = True)  # angle should be string

        by = (self.M*self.C*self.beta*self.gamma / self.Q) * (a * sp.pi / 180 / self.length)
        rho = self.M*self.C*self.beta*self.gamma / (self.Q * by)
        theta = l / rho
        C = sp.cos(theta)
        S = sp.sin(theta)

        M16 = rho * (1 - C) * (self.gamma / (self.gamma + 1))
        M26 = S * (self.gamma / (self.gamma + 1))
        M51 = self.f * (-S / (self.beta * self.C))
        M52 = self.f * (-rho * (1 - C) / (self.beta * self.C))
        M56 = self.f * (-rho * (l / rho - S) / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        mat = Matrix([[C, rho * S, 0, 0, 0, M16],
                      [-S / rho, C, 0, 0, 0, M26],
                      [0, 0, 1, l, 0, 0],
                      [0, 0, 0, 1, 0, 0],
                      [M51, M52, 0, 0, 1, M56],
                      [0, 0, 0, 0, 0, 1]])
      
        return mat

    def __str__(self):
        return f"Horizontal dipole magnet segment {self.length} m long (curvature) with an angle of {self.angle} degrees"

class dipole_wedge(lattice):
    def __init__(self, length, angle: float = 1, dipole_length: float = 0.0889, dipole_angle: float = 1.5, E0 = 0.51099, Q = 1.60217663e-19, M = 9.1093837e-31, E = 45, fringeType = 'decay'):
        super().__init__(length, E0, Q, M, E, fringeType)
        self.color = "lightgreen"
        self.angle = angle
        self.dipole_length = dipole_length
        self.dipole_angle = dipole_angle
        self.gap = 0.01  # hard-edge fringe field gap (m)
    
    def getSymbolicMatrice(self, numeric = False, length = None, angle = None):
        l = None
        a = None

        if length is None:
            l = self.length
        else:
            if numeric: l = length  # length should be number
            else: l = symbols(length, real = True)  # length should be string

        if angle is None:
            a = self.angle
        else:
            if numeric: a = angle  # angle should be number
            else: a = symbols(angle, real = True)  # angle should be string

        dipole_angle = self.dipole_angle
        dipole_length = self.dipole_length

        # Hard edge model for the wedge magnets
        By = (self.M*self.C*self.beta*self.gamma / self.Q) * (dipole_angle * sp.pi / 180 / dipole_length)
        R = self.M*self.C*self.beta*self.gamma / (self.Q * By)
        k = sp.Abs((self.Q * By / self.length) / (self.M * self.C * self.beta * self.gamma)) # Verify
        eta = (a * sp.pi / 180) * l / self.length # Verify
        E = (l * k) / ((sp.cos(eta)) ** 2)
        T = sp.tan(eta)
        M56 = self.f * (l / (self.C * self.beta * self.gamma * (self.gamma + 1)))

        mat = Matrix([[1, 0, 0, 0, 0, 0],
                      [-T / R, 1, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, -(T + E / 3) / R, 1, 0, 0],
                      [0, 0, 0, 0, 1, M56],
                      [0, 0, 0, 0, 0, 1]])
        
        return mat

    def __str__(self):
        return f"Horizontal wedge dipole magnet segment {self.length} m long (curvature) with an angle of {self.angle} degrees"
    


#NOTE: getSymbolicMatrice() must use all sympy methods and functions, NOT numpy
class beamline:
    class fringeField(lattice):
        B = 0 #  Teslas

        def __init__(self, length, fieldStrength,  E0=0.51099, Q=1.60217663e-19, M=9.1093837e-31, E=45):
            super().__init__(length, E0, Q, M, E)
            self.B = fieldStrength

        def __str__(self) -> str:
            return f"Fringe field segment {self.length} m long with a magnetic field of {self.B} teslas"
    
    def __init__(self, line = []):
        self.C = 299792458.0  # Speed of light in vacuum (m/s)
        self.q_e = 1.602176634e-19  # Elementary charge (C)
        self.m_e = 9.1093837139e-31  # Electron Mass (kg)
        self.m_p = 1.67262192595e-27  # Proton Mass (kg)
        self.m_amu = 1.66053906892E-27  # Atomic mass unit (kg)
        
        self.k_MeV = 1e-6 / self.q_e  # Conversion factor (MeV / J)

        #  [Mass (kg), charge (C), rest energy (MeV)]
        self.PARTICLES = {"electron": [self.m_e, self.q_e, (self.m_e * self.C ** 2) * self.k_MeV],
                          "proton": [self.m_p, self.q_e, (self.m_p * self.C ** 2) * self.k_MeV]}
        self.FRINGEDELTAZ = 0.01

        self.beamline = line
        self.totalLen = 0
        for seg in self.beamline:
            # seg.startPos = self.totalLen
            self.totalLen = self.totalLen + seg.length
            # seg.endPos = self.totalLen


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
            
     #  may use other interpolation methods (cubic, spline, etc)  
     # x = 0 = start of end of segment     
    def interpolateData(self, xData, yData, interval):
        rbf = interpolate.Rbf(xData, yData)
        totalLen = xData[-1] - xData[0]
        xNew = np.linspace(xData[0], xData[-1], math.ceil(totalLen/interval) + 1)
        yNew = rbf(xNew)
        return xNew, yNew
    '''
    ind: int
        The indice of the magnetic segment to create fringe
    '''
    # assumes the zlist doesnt start at 0, and that the first element is how far from the segment fringe measurement is
    def _addEnd(self, zList, magnetList, beamline, ind):
        
        driftLen = 0
        ind2 = ind
        while (ind2 != 0 and isinstance(beamline[ind2 - 1], driftLattice)):
            driftLen = driftLen + beamline[ind2 - 1].length
            ind2 -= 1
        
        i = 1
        fringeTotalLen = 0
        zList.insert(0,0)
        while (i < len(zList) and fringeTotalLen <= driftLen):
            fringeLen = zList[i] - zList[i-1]
            fringeTotalLen += fringeLen
            print(zList[i-1])

            if fringeTotalLen <= driftLen:
                fringeSeg = self.fringeField(fringeLen, magnetList[i-1])
                beamline.insert(ind, fringeSeg)
            
            i += 1
        
        while (fringeTotalLen > 0):
            if (beamline[ind-1].length <= fringeTotalLen):
                fringeTotalLen = fringeTotalLen - beamline[ind-1].length
                beamline.pop(ind-1)
                ind -= 1
            else:
                beamline[ind-1].length -= fringeTotalLen
                fringeTotalLen -= fringeTotalLen


            


    def createFringe(self, segment, fringeType, interval):
        pass
    
    # BEAMLINE OBJECT DOESNT CONTAIN THE BEAMLINE, ONLY TO PERFORM CALCULATIONS ON THE LINE

    # TEST: functiom ran on the class' beamline OBJECT, NOT beamline LIST

    #  Fringe fields can only exist overlapping drift segments for now
    #  assume all fringe fields overlap drift segments
    def reconfigureLine(self, interval = None):
        if interval is None:
            interval = self.FRINGEDELTAZ

        beamline = self.beamline

        # zPos = []
        # magStrength = []
        # i = 0
        # while (i <= self.totalLen):
        #     zPos.append(i)
        #     magStrength.append(0)
        #     i = i + interval
        
        for segment in beamline:
            if isinstance(segment.fringeType, list):
                xData = segment.fringeType[0]
                yData = segment.fringeType[1]
                xNew, yNew = self.interpolateData(xData, yData, interval)

                self._addEnd(yNew)

                # segmentZ = []
                # for i in zPos: segmentZ.append(i)
