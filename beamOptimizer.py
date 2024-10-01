#Author: Christian Komo

'''
helpful resources
https://www.youtube.com/watch?v=G0yP_TM-oag
https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
https://docs.scipy.org/doc/scipy/reference/optimize.html
gpt"can you xplain the parameters of a constraint paramter in scipy.minimize"
'''
#  NOTE: nelder-mead method doesn't work if starting search point is zero (if variablesValues = [0,0,0...])

import scipy.optimize as spo
from beamline import *
from ebeam import beam
import numpy as np
from schematic import *
import timeit
import matplotlib.pyplot as plt

class beamOptimizer():
    def __init__(self, beamline, segmentVar: dict, 
                 stdxend, stdyend, method, 
                 matrixVariables, startPoint = {}, xWeight = 1, yWeight = 1, objectives = {}):
        '''
        Constructor for the optimizer object. Object is used to optimize the electric current values
        for quadruples in an accelerator beamline in order that desired particle x and y positional spread may be
        acheived

        Parameters
        ----------
        beamline: list[beamline]
            list of beamline objects representing accelerator beam

        segmentVar: list(int)
            list of segmentVar representing the interval in the beamline we want to optimize 

        stdxend: float
            final standard deviation of particles' x position that we are targeting
        stdyend: float
            final standard deviation of particles' y position that we are targeting
        start: list[float]
            list of electrical current values for minimizing function to start optimizing from. 
            Each nth value in list HAS corresponds to the nth qpd/qfd.
        method: str
            name of optimization method to use for desired values
        matrixVariables: np.array(list[float][float])
            2D numPy list of particle elements
        xWeight: float, optional
            number giving the algorithm more or less bias towards minimizing standard deviation goal difference for x position.
            weight > 1 means more bias, weight < 1 means less bias
        yWeight: float, optional
            number giving the algorithm more or less bias towards minimizing standard deviation goal difference for y position.
            weight > 1 means more bias, weight < 1 means less bias

        '''
        self.stdxend = stdxend
        self.stdyend = stdyend
        self.matrixVariables = matrixVariables
        self.beamline = beamline
        self.method = method
        self.xWeight = xWeight
        self.yWeight = yWeight
        self.objectives = objectives


        self.plotChiSquared = []
        self.plotIterate = []
        self.iterationTrack = 0

        self.segmentVar = segmentVar
        checkSet = set()
        self.variablesToOptimize = []
        for item in segmentVar:
            if (item < 0 or item >= len(beamline)):
                raise IndexError
            varItem = segmentVar.get(item)[0]
            if varItem not in checkSet:
                self.variablesToOptimize.append(varItem)
                checkSet.add(varItem)

        self.variablesValues = [] 
        self.bounds = []
        self.trackVariables = []
        for i in self.variablesToOptimize:
            self.variablesValues.append(1) 
            self.bounds.append((None, None))
            self.trackVariables.append([])
        
        for var in startPoint:
            index = self.variablesToOptimize.index(var)
            if "start" in startPoint.get(var): self.variablesValues[index] = startPoint.get(var).get("start")
            if "bounds" in startPoint.get(var): self.bounds[index] = startPoint.get(var).get("bounds")

    def getXStd(particles):
        return np.std(particles[:,0])
        
    def getYStd(particles):
        return np.std(particles[:,2])
        

    def _optiSpeed(self, variableVals):
        '''
        Simulates particle movement through a beamline, calculates positional standard deviation, 
        and returns a chi-squared accuracy statistic. Function to call on for standard deviation
        optimization.

        Parameters
        ----------
        current: list[float]
            test current value(s) for quadruples (IT IS ASSUMED THAT THE EACH nth 
            ELEMENT OF current CORRESPONDS TO THE nth NUMBER OF a QPF/QPD)

        Returns
        -------
        difference: float
            chi-squared statistic, it is the combined squared difference between 
            target standard deviation and actual standard deviation devided by the 
            target standard deviation. Weighted bias for x vs y stats also accounted for.
        '''
        particles = self.matrixVariables
        segments = self.beamline


        for i in range(len(segments)):
            if i in self.segmentVar:
                yFunc = self.segmentVar.get(i)[1]
                varIndex = self.variablesToOptimize.index(self.segmentVar.get(i)[0]) #  Get the index of the variable to use with 
                objCurrent = yFunc(variableVals[varIndex])
                param = self.segmentVar.get(i)[2]
                particles = np.array(segments[i].useMatrice(particles, **{param: objCurrent}))

                if i in self.objectives:
                    chiPiece = 
            else:
                particles = np.array(segments[i].useMatrice(particles))  
                
        stdx = np.std(particles[:,0])
        stdy = np.std(particles[:,2])


        difference = np.sqrt(((stdx-self.stdxend)**2)*self.xWeight*(1/self.stdxend) + ((stdy-self.stdyend)**2)*self.yWeight*(1/self.stdyend))
        self.plotChiSquared.append(difference)
        self.plotIterate.append((self.iterationTrack) + 1)
        self.iterationTrack = self.iterationTrack + 1
        self.trackVariables[0].append(stdx)
        self.trackVariables[1].append(stdy)
        # print(difference)  #for testing
        return difference
    
    def calc(self, plot = False):
        '''
        optimizes beamline quadruple current values so particles' positional standard
        deviation is as close to target as possible

        Returns
        -------
        result: OptimizeResult
            Object containing resulting information about optimization process and results
        '''

        # result = spo.minimize(self._optiSpeed, self.start, options={"disp": True}, method = self.method)
        result = spo.minimize(self._optiSpeed, self.variablesValues, method = self.method, bounds=self.bounds, options={'disp':True})

        if plot:
            fig, ax1 = plt.subplots()

            chiLine, =ax1.plot(self.plotIterate, self.plotChiSquared, label = 'Chi Squared', color = 'green')
            ax1.set_yscale('log')
            ax1.set_ylabel('Chi-Squared Difference', color = 'green')
            ax1.tick_params(axis='y', labelcolor = 'green')
            ax1.spines['left'].set_color('green')
            ax1.spines['left'].set_linewidth(1)

            ax2 = ax1.twinx()
            xStdLine, = ax2.plot(self.plotIterate, self.trackVariables[0], label = 'xstd', color = 'red')
            yStdLine, =ax2.plot(self.plotIterate, self.trackVariables[1], label = 'ystd', color = 'blue')
            ax2.set_ylabel('x and y standard deviation')
            

            plt.legend(handles = [chiLine, xStdLine, yStdLine], loc = 'upper right')
            plt.show()

        #  WIP: add results in a clean format with values of each segment
        #  Potentially add this to a seperate tinker window?
        for key in self.segmentVar:
            variable = self.segmentVar.get(key)[0]
            index = self.variablesToOptimize.index(variable)
            yFunc = self.segmentVar.get(key)[1]
            newVal = yFunc(result.x[index])
            self.segmentVar.get(key).append(newVal)


        return result
   
    def testSpeed(self, iterations):
        '''
        Test the speed of an optimization algorithm (more iterations = more accurate)

        Parameters
        ----------
        iterations: float
            Number of times to run calc() function with 

        Returns
        -------
        timeResult: float
            average number of seconds to execute calc() function once
        '''
        timeResult = (timeit.timeit(self.calc,number = iterations))/iterations
        return timeResult
    
    def testFuncEval(self, iterations):
        '''
        Test the number of function evalutions an optimization algorithm performs (more iterations = more accurate)

        Parameters
        ----------
        iterations: float
            Number of times to run calc() function with 

        Returns
        -------
        timeResult: float
            average number of function evalutions per calc() call
        evalType: str
            type of evaluation that is being returned
        '''
        iterationsTotal = 0
        evalType = ''
        for i in range(iterations):
            result = self.calc()
            try: 
                iterationsTotal += result.nfev
                evalType = "# of function evaluations"
            except AttributeError: 
                try: 
                    iterationsTotal += result.njev
                    evalType = "# of Jacobian evaluations"
                except AttributeError: 
                    iterationsTotal += result.nhev
                    evalType = "# of Hessian evaluations"
        return iterationsTotal/iterations, evalType   
        
    def testFuncIt(self, iterations):
        '''
        Test the number of iterations of an optimization algorithm (more iterations = more accurate)

        Parameters
        ----------
        iterations: float
            Number of times to run calc() function with 
        
        Returns
        -------
        avIterate: float
            average number of function iterations per calc() call
        '''
        iterationsTotal = 0
        for i in range(iterations):
            result = self.calc()
            try: iterationsTotal += result.nit
            except AttributeError: 
                raise AttributeError("algorithm does not track iterations")
        avIterate = iterationsTotal/iterations   
        return avIterate
    
