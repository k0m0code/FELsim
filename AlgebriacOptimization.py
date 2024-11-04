from beamline import *
from ebeam import *
import sympy as sp

class AlgebriacOpti():
    def __init__(self):
        self.DDOF = 1
        self.m = None
        self.sigmai = None
        pass

    def getSigmai(self, particles):
        ebeam = beam()
        dist_avg, dist_cov, twiss = ebeam.cal_twiss(particles, ddof=self.DDOF)
        sigmaI = [[0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0]]
        print(twiss)
        for i, axis in enumerate(['x', 'y', 'z']):
            ax = twiss.loc[axis]
            epsilon = ax["$\epsilon$ ($\pi$.mm.mrad)"]
            alphaEpsilon = ax[r"$\alpha$"]*epsilon
            beta = ax[r"$\beta$ (m)"]
            gamma = ax[r"$\gamma$ (rad/m)"]
    
            sigmaI[i*2][i*2] = beta*epsilon
            sigmaI[i*2][i*2+1] = -(alphaEpsilon)
            sigmaI[i*2+1][i*2] = -(alphaEpsilon)
            sigmaI[i*2+1][i*2+1] = gamma*epsilon
        return sp.Matrix(sigmaI)

         


    def getM(self, beamline: list, xVar: dict):
        resultArr = None
        try:
             resultArr = beamline[-1].getSymbolicMatrice(**xVar[len(beamline) - 1])
        except KeyError:
            resultArr = beamline[-1].getSymbolicMatrice()
        i = len(beamline) - 2
        while (i >= 0):
                try:
                    resultArr = resultArr*beamline[i].getSymbolicMatrice(**xVar[i])
                except KeyError:
                    resultArr = resultArr*beamline[i].getSymbolicMatrice()
                i = i - 1
        return resultArr
    
    def getSigmaF(self, m, sigmaI):
         mTransposed = m.T
         return m*sigmaI*mTransposed
         
    


                


# import numpy as np
# from sympy import symbols, Matrix

# # Define your symbolic variables
# a, b, c, d = symbols('a b c d')

# # Create symbolic matrices
# A = Matrix([[a, b], 
#             [c, d]])
# B = Matrix([[1, 2],
#             [3, 4]])

# # Perform matrix multiplication
# C = A * B

# print(C)  # This will print the resulting symbolic matrix
