#Author: Christian Komo

'''
helpful resources
https://www.youtube.com/watch?v=G0yP_TM-oag
https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
https://docs.scipy.org/doc/scipy/reference/optimize.html
gpt"can you xplain the parameters of a constraint paramter in scipy.minimize"
'''
#  NOTE: nelder-mead method doesn't work if starting search point is zero (if variablesValues = [0,0,0...]
#  NOTE: After testing, COBYLA and some other methods may try a test value of 0 for current, length, etc, that may throw back a difference of NAN
#        (because beamline object divides by 0). Have to figure out how to bound variable values automatically so computer doesn't use weird values like 0, 
#        negative numbers (YOU MUST USE bounds for now so program doesn't use negative/zero numbers

import scipy.optimize as spo
from beamline import *
from ebeam import beam
import numpy as np
from schematic import *
import timeit
import matplotlib.pyplot as plt

# NOTE: EACH BEAMOPTIMIZER OBJECT AFTER INSTANTIATION SHOULD ONLY BE USED TO 
# RUN A CALC() FUNCTION ONE TIME, CODE HAS NOT BEEN MODIFIED YET BEYOND ONE CALC() FUNCTION CALL

# For each indice of the beam segment in parameter, no proper error handling yet for invalid values, repeating indices with same objective, have to test...

#  Currently can only optimize one variable for each segment indice, have to implement more than one variable inthe future

# TODO: some objectives would like value goal to be zero, but we cant divide by zero in our chi difference calc. Find a way



class beamOptimizer():
    def __init__(self, beamline, segmentVar: dict, method, 
                 matrixVariables,objectives, startPoint = {}):
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
        self.matrixVariables = matrixVariables
        self.beamline = beamline
        self.method = method
        
        ebeam = beam()
        self.OBJECTIVEMETHODS = {"std": ebeam.std,"epsilon":ebeam.epsilon,"alpha":ebeam.alpha,"beta":ebeam.beta,"gamma":ebeam.gamma,"phi":ebeam.phi} #  Methods included in class to return staistical information

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

        
        self.plotChiSquared = []
        self.plotIterate = []
        self.iterationTrack = 0
        self.trackVariables = {}

        #  NOTE: "measure" has to be a function call that returns a single value with a parameter of a 2d list of particles, and each indice can only appear once as a key
        self.objectives = objectives
        for key, value in self.objectives.items():
            for goal in value:
                if goal["measure"][1] in self.OBJECTIVEMETHODS:
                    goal["measure"][1] = self.OBJECTIVEMETHODS[goal["measure"][1]]
                self.trackVariables.update({"indice " + str(key) + ": " + goal["measure"][0] + " "  + goal["measure"][1].__name__: []})


        self.variablesValues = [] 
        self.bounds = []
        for i in self.variablesToOptimize:
            self.variablesValues.append(1) 
            self.bounds.append((None, None))
        for var in startPoint:
            index = self.variablesToOptimize.index(var)
            if "start" in startPoint.get(var): self.variablesValues[index] = startPoint.get(var).get("start")
            if "bounds" in startPoint.get(var): self.bounds[index] = startPoint.get(var).get("bounds")
    
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
        chiPieces = []


        for i in range(len(segments)):
            if i in self.segmentVar:
                yFunc = self.segmentVar.get(i)[1]
                varIndex = self.variablesToOptimize.index(self.segmentVar.get(i)[0]) #  Get the index of the variable to use with 
                newValue = yFunc(variableVals[varIndex])
                param = self.segmentVar.get(i)[2]
                particles = np.array(segments[i].useMatrice(particles, **{param: newValue}))
            else:
                particles = np.array(segments[i].useMatrice(particles))  
            if i in self.objectives:
                for goalDict in self.objectives[i]:
                    stat = (goalDict["measure"][1](particles, goalDict["measure"][0]))
                    chiPieces.append((((stat-goalDict["goal"])**2)*goalDict["weight"])/goalDict["goal"])
                    stringForm = "indice " + str(i) + ": " + goalDict["measure"][0] + " " + goalDict["measure"][1].__name__
                    self.trackVariables[stringForm].append(stat)

        difference = 0
        difference = np.sqrt(np.sum(chiPieces))

        self.plotChiSquared.append(difference)
        self.plotIterate.append((self.iterationTrack) + 1)
        self.iterationTrack = self.iterationTrack + 1

        # print("diff:" + str(difference))  #for testing
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
        result = spo.minimize(self._optiSpeed, self.variablesValues, method = self.method, bounds=self.bounds, options={'disp':True})
        #  Alpha parameters look a little weird when plotting their minimization
        if plot:
            fig, ax1 = plt.subplots()
            handles = []

            chiLine, =ax1.plot(self.plotIterate, self.plotChiSquared, label = 'Chi Squared', color = 'green')
            ax1.set_yscale('log')
            ax1.set_ylabel('Chi-Squared Difference', color = 'green')
            ax1.tick_params(axis='y', labelcolor = 'green')
            ax1.spines['left'].set_color('green')
            ax1.spines['left'].set_linewidth(1)
            handles.append(chiLine)

            ax2 = ax1.twinx()
            mini = 0
            for i, key in enumerate(self.trackVariables):
                valAx, = ax2.plot(self.plotIterate, self.trackVariables[key], label = key)
                handles.append(valAx)
                tempMin = abs(min(self.trackVariables[key]))
                if i == 0 or mini>tempMin:
                    mini = tempMin
            print(str(10**(np.ceil(np.log10(mini)))))
            print("min" + str(mini))
            ax2.set_yscale('symlog', linthresh=10**(np.ceil(np.log10(mini))))
            ax2.set_ylabel('Objective functions')
            

            plt.legend(handles = handles, loc = 'upper right')
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
    
