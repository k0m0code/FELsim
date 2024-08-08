#Author: Christian Komo


'''
https://www.youtube.com/watch?v=G0yP_TM-oag
https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
https://docs.scipy.org/doc/scipy/reference/optimize.html
gpt"can you xplain the parameters of a constraint paramter in scipy.minimize"
'''
import scipy.optimize as spo
from beamline import *
import numpy as np



class beamOptimizer():
    def __init__(self, matrixVariables, beamline, indices: tuple, stdxend, stdyend):
        self.stdxend = stdxend
        self.stdyend = stdyend
        self.matrixVariables = matrixVariables
        self.beamline = beamline
        self.indices = indices

    def updateIndices(self, indices: tuple):
        self.indices = indices

    def func(self, current):
        segment = qpfLattice(current = current,length = 80)
        endValues = np.array(segment.useMatrice(values = self.matrixVariables))
        stdx = np.std(endValues[:,0])
        stdy = np.std(endValues[:,2])
        return np.sqrt((stdx-self.stdxend)**2+(stdy-self.stdyend)**2)
    
    def func2(self, current):
        particles = self.matrixVariables
        segments = self.beamline[self.indices[0]:self.indices[1]]
        ii = 0
        for i in range(len(segments)):
            if isinstance(segments[i], qpdLattice) or isinstance(segments[i], qpfLattice):
                particles = np.array(segments[i].useMatrice(particles, current = current[ii]))
                ii += 1
            else:
                particles = np.array(segments[i].useMatrice(particles))
        stdx = np.std(particles[:,0])
        stdy = np.std(particles[:,2])
        return np.sqrt((stdx-self.stdxend)**2+(stdy-self.stdyend)**2)
    
    def calc(self, startx):
        # constrain = 
        result = spo.minimize(self.func2, startx, options={"disp": True})
        return result