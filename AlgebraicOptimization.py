from beamline import *
from ebeam import *
import sympy as sp

class AlgebraicOpti():
    def __init__(self):
        self.DDOF = 1
        self.m = None
        self.sigmai = None

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
    
    def getSigmaF(self, m, sigmaI):
         mTransposed = m.T
         return m*sigmaI*mTransposed
        
    def findObj(self, beamline, xVal, objec, startParticles = None):
        sigi = None
        if not startParticles is None:
              sigi = self.getDistSigmai(startParticles)
        else:
            objList = []
            for i, axis in enumerate(['x','y','z']):
                objList.append(objec[axis])
            sigi = self.getTwissSigmai(objList[0],objList[1],objList[2])
        mMat = self.getM(beamline, xVal)
        sigF = self.getSigmaF(mMat,sigi)
        for i in range(len(sigF)):
             sigF[i] = sigF[i] - sigi[i]
        return sigF # return the objective functions, solutions at zeros
         

    
