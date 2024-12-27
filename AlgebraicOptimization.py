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
        self.MINCHECK = 0.00001
        self.CURRENTRANGE = 10
        self.ROOTINTERVAL = 0.5
        self.DOMAIN = (0.00001,10)
        

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
    
    def getSigmaF(self, transferM, sigmaI):
        '''
        returns linear equations representing final beam matrix
        '''  
        mTransposed = transferM.T
        return transferM*sigmaI*mTransposed
    

    #  LINEAR EQUATIONs TO OPTIMIZE TO STARTING y TWISS CONDITIONS AS POSSIBLE 
    #  assumed that particles/ twiss are begiging beamline parameters
    def findSymmetricObjective(self, beamline, xVar, startParticles = None, twiss = None, latex = False):
        '''
        returns the transformed beam matrix from beamline and twiss/particle parameters.

        Parameters
        ----------
        beamline: list[beamline]
            list of beamline objects representing accelerator beam
        xVar: dict
            dictionary of beamline element indices and their parameter values to optimize
        startParticles = np.array(list[float][float])
            2D numPy list of particle elements
        twiss: dict
            twiss values for each dimensional plane
        latex: bool
            whether to return equations in Latex form

        Returns
        -------
        latexList: list[list[str]]
            6x6 2D list of equations in latesx form represented as strings
        sigObg: sympy.Matrix
            6x6 2D Symypy matrix of equations representing transformed beam matrix
        '''
        sigi = None
        if not startParticles is None:
              sigi = self.getDistSigmai(startParticles)
        elif not twiss is None:
            objList = []
            for ind, axis in enumerate(['x','y','z']):
                objList.append(twiss[axis])
            sigi = self.getTwissSigmai(objList[0],objList[1],objList[2])
        else:
             raise ValueError("Please enter twiss or startParticles parameter")
        if (not twiss is None) and (not startParticles is None):
             raise ValueError("Please enter one parameter for either twiss or startParticles only")
        mMat = self.getM(beamline, xVar)
        sigObg = self.getSigmaF(mMat,sigi)
        for i in range(len(sigObg)):
             sigObg[i] = sigObg[i] - sigi[i]

        if latex:
            latexList = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
            for i in range(len(latexList)):
                for ii in range(len(latexList[i])):
                    latexList[i][ii] = sp.latex(sigObg[i,ii])
            return latexList
        return sigObg  # return the objective functions, solutions at zeros
         
    #NOTE: only created to find roots for only one variable that exists in the 
    #equation
    def getRootsUni(self, equation):
        '''
        returns the roots of a complex equation using the multi-start method,
        function to be used when finding segment parameter values for a sigmaf equation. 
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
        rootSet = set()
        ind = self.DOMAIN[0]
        tempSet = equation.free_symbols
        var = tempSet.pop()
        total_intervals = self.DOMAIN[1]
        with tqdm(total=total_intervals, desc="Finding roots...",
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]") as pbar:
            while ind <= self.DOMAIN[1]:
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
    
    # create class variable to have a obvious x-y = 0 line 
    def getRootsMulti(self, equation):
        '''
        returns the roots of a complex implicit equation using the multi-start method, that is, where
        the implicit equation equals x - y = 0.
        Function to be used when finding segment parameter values for a sigmaf equation. 
        Function intended for bivariate implicit equations
        
        Parameters
        ----------
        equation: sympy.add
            sympy equation to find roots at

        Returns
        -------
        rootList: list[tuple]
            estimated list of zero pairs of equation in specified interval
        '''
        ans = None
        rootSet = set()
        sett = equation.free_symbols
        x = sett.pop()
        y = sett.pop()
        checkLine = sp.Eq(x-y, 0)
        ind = self.DOMAIN[0]
        total_intervals = self.DOMAIN[1]
        with tqdm(total=total_intervals, desc="Finding roots...",
                  bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]") as pbar:
            while ind <= self.DOMAIN[1]:
                try:
                    ans = sp.nsolve((equation, checkLine), (x, y), (ind, ind))
                    if (ans[0] > 0 and ans[1] > 0):
                        tup = tuple(tuple(row) for row in ans.tolist())
                        rootSet.add(tup)
                except ValueError:
                    pass
                pbar.update(self.ROOTINTERVAL)
                ind = ind + self.ROOTINTERVAL
        rootList = list(rootSet)
        return rootList
