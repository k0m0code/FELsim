from beamOptimizer import beamOptimizer
from ebeam import beam
from schematic import draw_beamline
from beamline import *
import numpy as np

ebeam = beam()
schem = draw_beamline()

I = 2.28
sec1 = driftLattice(0.5)
sec2 = qpfLattice(current = I)
sec3 = driftLattice(0.25)
sec4 = qpdLattice(current = I)
sec5 = driftLattice(0.25)
sec6 = qpfLattice(current = I)
sec7 = driftLattice(0.25)
sec8 = qpdLattice(current = I)
sec9 = driftLattice(0.50)
sec10 = dipole(length=0.0889, angle=1.5)

line = [sec1,sec2,sec3,sec4,sec5,sec6,sec7,sec8,sec9, sec10]

beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,1,1],1000)


vals = {0: ["A", lambda num: num, "length"],
        1: ["I", lambda num:num, "current"],
        3: ["I", lambda num:num, "current"],
        5: ["I", lambda num:num, "current"],
        7: ["I", lambda num:num, "current"],}

starting = {"I": {"bounds": (0,10), "start": 1}, 
            "A": {"start": 0.5}}

def alpha(particles):
    return beam.getXYZ(particles)

objectives = {9: {"measure": "xStd", "goal": 5, "weight": 9}, 5: {"measure": "yStd", "goal": 5, "weight": 9}}

valTest =       {0: {"variable": "A", "relationship": lambda num: num, "optimized": "length"},
                1: ["I", lambda num:num, "current"],
                3: ["I", lambda num:num, "current"],
                5: ["I", lambda num:num, "current"],
                7: ["I", lambda num:num, "current"],}

matrixVariables = ebeam.gen_6d_gaussian(0,[1,.2,1,0.2,1,1],1000)
test = beamOptimizer(line, vals, "COBYLA", matrixVariables, startPoint= starting, objectives=objectives)
print(test.trackVariables)





# test = beamOptimizer(line, (1,3,5),1,1, [A,B], "Nelder-Mead", matrixVariables)


result = test.calc(plot=True)


# print("speed " + str(test.testSpeed(10)))
# evals, evalType = test.testFuncEval(10)
# print("evals " + str(evals) + ", " + evalType)
# print("iterations " + str(test.testFuncIt(10)))

beam_dist = test.matrixVariables
schem = draw_beamline()
line[0].length = result.x[0]
line[1].current = result.x[1]
line[3].current = result.x[1]
line[5].current = result.x[1]
line[7].current = result.x[1]
schem.plotBeamPositionTransform(beam_dist, line, 0.05)

print("Current" + str(result.x))
print("Chi Squared:" + str(result.fun))
print("xstd; " + str(np.std(schem.matrixVariables[:, 0])))
print("ystd: " + str(np.std(schem.matrixVariables[:, 2])))
