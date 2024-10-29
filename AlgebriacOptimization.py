from beamline import *
from ebeam import *

class AlgebriacOpti():
    def __init__(self):
        self.DDOF = 1
        self.m = None
        self.sigmai = None
        pass

    def getSigmai(self, particles):
         ebeam = beam()
         dist_avg, dist_cov, twiss = ebeam.cal_twiss(particles, ddof=self.DDOF)

    def getM(self, beamline):
            resultArr = beamline[-1]
            i = len(beamline) - 1
            # while (i >= 0):
            #     #  resultArr = 


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


