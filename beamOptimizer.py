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
import numpy as np
from schematic import *
import timeit

class beamOptimizer():
    def __init__(self, beamline, optimVariables: dict, 
                 stdxend, stdyend, method, 
                 matrixVariables, startPoint = {}, xWeight = 1, yWeight = 1):
        '''
        Constructor for the optimizer object. Object is used to optimize the electric current values
        for quadruples in an accelerator beamline in order that desired particle x and y positional spread may be
        acheived

        Parameters
        ----------
        beamline: list[beamline]
            list of beamline objects representing accelerator beam

        optimVariables: tuple(int)
            tuple of optimVariables representing the interval in the beamline we want to optimize 

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

        self.optimVariables = optimVariables
        checkSet = set()
        self.variablesOptimize = []
        for item in optimVariables:
            varItem = optimVariables.get(item)[0]
            if varItem not in checkSet:
                self.variablesOptimize.append(varItem)
                checkSet.add(varItem)

        self.variablesValues = [1 for i in self.variablesOptimize] 
        self.bounds = [(None, None) for i in self.variablesOptimize]
        for var in startPoint:
            index = self.variablesOptimize.index(var)
            if "start" in startPoint.get(var): self.variablesValues[index] = startPoint.get(var).get("start")
            if "bounds" in startPoint.get(var): self.bounds[index] = startPoint.get(var).get("bounds")

        print(self.variablesOptimize)
        
        self.method = method
        self.xWeight = xWeight
        self.yWeight = yWeight

    def _optiSpeed(self, current):
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
            if i in self.optimVariables:
                yFunc = self.optimVariables.get(i)[1]
                varIndex = self.variablesOptimize.index(self.optimVariables.get(i)[0]) #  Get the index of the variable to use with 
                objCurrent = yFunc(current[varIndex])
                param = self.optimVariables.get(i)[2]
                particles = np.array(segments[i].useMatrice(particles, **{param: objCurrent}))
            else:
                particles = np.array(segments[i].useMatrice(particles))  
                
        stdx = np.std(particles[:,0])
        stdy = np.std(particles[:,2])
        difference = np.sqrt(((stdx-self.stdxend)**2)*self.xWeight*(1/self.stdxend) + ((stdy-self.stdyend)**2)*self.yWeight*(1/self.stdyend))
        # print(difference)  #for testing
        return difference
    








    













    
    def calc(self):
        '''
        optimizes beamline quadruple current values so particles' positional standard
        deviation is as close to target as possible

        Returns
        -------
        result: OptimizeResult
            Object containing resulting information about optimization process and results
        '''

        # result = spo.minimize(self._optiSpeed, self.start, options={"disp": True}, method = self.method)
        result = spo.minimize(self._optiSpeed, self.variablesValues, method = self.method, bounds=self.bounds)
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
    
