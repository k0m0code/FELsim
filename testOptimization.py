from beamOptimizer import beamOptimizer
from ebeam import beam
from schematic import draw_beamline
from beamline import *
import numpy as np

ebeam = beam()

sec1 = driftLattice(200)
sec2 = qpfLattice(current = 32)
sec3 = driftLattice(100)
sec4 = qpdLattice(current = 32)
line = [sec1,sec2, sec3, sec4]

test = beamOptimizer(line, (0,4),50,90, [60,40], method = "COBYLA")
result = test.calc()
# print("speed " + str(test.testSpeed(10)))
# evals, evalType = test.testFuncEval(10)
# print("evals " + str(evals) + ", " + evalType)
# print("iterations " + str(test.testFuncIt(10)))

beam_dist = test.matrixVariables
schem = draw_beamline()
line[1].current = result.x[0]
line[3].current = result.x[1]
schem.plotBeamPositionTransform(beam_dist, line, 10)
print("Current" + str(result.x))
print("Chi Squared:" + str(result.fun))
print("xstd; " + str(np.std(schem.matrixVariables[:, 0])))
print("ystd: " + str(np.std(schem.matrixVariables[:, 2])))
