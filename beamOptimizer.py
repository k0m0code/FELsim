#Author: Christian Komo

'''
helpful resources
https://www.youtube.com/watch?v=G0yP_TM-oag
https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
https://docs.scipy.org/doc/scipy/reference/optimize.html
gpt"can you xplain the parameters of a constraint paramter in scipy.minimize"
'''
import scipy.optimize as spo
from beamline import *
import numpy as np
from ebeam import beam
from schematic import *
import csv
import timeit

class beamOptimizer():
    def __init__(self, beamline, indices: tuple, 
                 stdxend, stdyend, start, method):
        '''
        Constructor for the optimizer object. Object is used to optimize the electric current values
        for quadruples in an accelerator beamline in order that desired particle x and y positional spread may be
        acheived

        Parameters
        ----------
        matrixvairables: np.array(list[float][float])
            2D numPy list of particle elements
        beamline: list[beamline]
            list of beamline objects representing accelerator beam
        indices: tuple(int)
            tuple of indices representing the interval in the beamline we want to optimize 
        stdxend: float
            final standard deviation of particles' x position that we are targeting
        stdyend: float
            final standard deviation of particles' y position that we are targeting
        start: list[float]
            list of electrical current values for minimizing function to start optimizing from. 
            Each nth value in list HAS corresponds to the nth qpd/qfd.
        method: str
            name of optimization method to use for desired values

        '''
        self.stdxend = stdxend
        self.stdyend = stdyend
        self.matrixVariables = None
        self.beamline = beamline
        self.indices = indices
        self.start = start
        self.method = method
    
    def _optiSpeed(self, current):
        '''
        simulates particle movement through a beamline, calculates positional standard deviation, 
        and returns an accuracy statistic

        Parameters
        ----------
        current: list[float]
            test current value(s) for quadruples (IT IS ASSUMED THAT THE EACH ELEMENT OF current CORRESPONDS TO THE NTH NUMBER OF QPF/QPD)

        returns
        -------
        difference: float
            the combined squared difference between target standard deviaiton and actual standard deviation
        '''
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
        #  difference = (((stdx-self.stdxend)**2)/self.stdxend) + (((stdy-self.stdyend)**2)/self.stdyend)
        difference = (stdx-self.stdxend)**2 + (stdy-self.stdyend)**2
        return difference
    
    def calc(self):
        '''
        Generates random particles and optimizes beamline variables so particles' positional standard
        deviation is as close to target as possible

        Returns
        -------
        result: OptimizeResult
            Object containing resulting information about optimization process and results
        '''
        ebeam = beam()
        self.matrixVariables = ebeam.gen_6d_gaussian(0,[1,.1,1,0.1,1,1],1000)

        # result = spo.minimize(self._optiSpeed, self.start, options={"disp": True}, method = self.method)
        result = spo.minimize(self._optiSpeed, self.start, method = self.method)
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
            type of evaluations that is being returned
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
