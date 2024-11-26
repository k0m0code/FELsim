from beamline import *
from ebeam import *
import sympy as sp
import sys
from tqdm import tqdm

'''
helpful resources
https://stackoverflow.com/questions/38104025/sympy-solveset-returns-conditionset-when-it-should-be-returning-a-numerical
'''

class AlgebraicOpti():
    def __init__(self):
        self.DDOF = 1
        self.CURRENTRANGE = 10
        self.ROOTINTERVAL = 0.5
        self.MINCHECK = 0.00001

    def getDistSigmai(self, particles):
        ebeam = beam()
        dist_avg, dist_cov, twiss = ebeam.cal_twiss(particles, ddof=self.DDOF)
        sigmaI = [[0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0]]
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
    
    #[epsilon, alpha, beta, gamma] for each axis
    def getTwissSigmai(self, xTwiss, yTwiss, zTwiss):
        twiss = [xTwiss,yTwiss,zTwiss]
        sigmaI = [[0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,0]]
        for i in range(3):
             ax = twiss[i]
             epsilon = ax[0]
             alphaEpsilon = -(ax[1]*epsilon)
             betaEpsilon = ax[2]*epsilon
             gammaEpsilon = ax[3]*epsilon
             sigmaI[i*2][i*2] = betaEpsilon
             sigmaI[i*2][i*2+1] = -(alphaEpsilon)
             sigmaI[i*2+1][i*2] = -(alphaEpsilon)
             sigmaI[i*2+1][i*2+1] = gammaEpsilon
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
    
    '''
    returns linear equations representing propogation of the beam
    '''  
    def getSigmaF(self, m, sigmaI):
         mTransposed = m.T
         return m*sigmaI*mTransposed
    

    #  LINEAR EQUATIONs TO OPTIMIZE TO STARTING y TWISS CONDITIONS AS POSSIBLE 
    
    def findObj(self, beamline, xVal, objec = None, startParticles = None):
        sigi = None
        if not startParticles is None:
              sigi = self.getDistSigmai(startParticles)
        elif not objec is None:
            objList = []
            for ind, axis in enumerate(['x','y','z']):
                objList.append(objec[axis])
            sigi = self.getTwissSigmai(objList[0],objList[1],objList[2])
        else:
             raise ValueError("Please enter objec or startParticles parameter")
        if (not objec is None) and (not startParticles is None):
             raise ValueError("Please enter one parameter for either objec or startParticles only")
        mMat = self.getM(beamline, xVal)
        sigObg = self.getSigmaF(mMat,sigi)
        for i in range(len(sigObg)):
             sigObg[i] = sigObg[i] - sigi[i]


        return sigObg # return the objective functions, solutions at zeros
         
    #NOTE: only created to find roots for only one variable that exists in the 
    #equation
    '''
    returns the roots of a complex equation using the multi-start method,
    function to be used when finding current amplitude value for a sigmaf equation. 
    Function intended for univariate equations
    
    Parameters
    ----------
    equation: sympy.add
        sympy equation to find roots at

    Returns
    -------
    rootList: list[float]
        estimated list of zeros of equation in specified interval
    '''
    def getRootsUni(self, equation):
            rootSet = set()
            ind = self.MINCHECK
            tempSet = equation.free_symbols
            var = tempSet.pop()
            total_intervals = self.CURRENTRANGE
            with tqdm(total=total_intervals, desc="Finding roots...",
                  bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]") as pbar:
                while ind <= self.CURRENTRANGE:
                    try:
                        val = sp.nsolve(equation, var ,ind)
                        if val > 0:
                            rootSet.add(val)
                    except ValueError:
                        pass
                    pbar.update(self.ROOTINTERVAL)
                    ind = ind + self.ROOTINTERVAL
            rootList = list(rootSet)
            return rootList
    
    # def getRootsMulti(self, equation):
         
