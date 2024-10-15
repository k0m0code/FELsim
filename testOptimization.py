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

beamtype = beamline()
pBeam = beamtype.changeBeamType(line, "electron", 55)

beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,1,1],1000)
schem.plotBeamPositionTransform(beam_dist, line, 0.05)


vals = {1: ["B", lambda num:num, "current"],
        3: ["I", lambda num:num, "current"],
        5: ["I", lambda num:num, "current"],
        7: ["A", lambda num:num, "current"]}

starting = {"I": {"bounds": (0.00001,10), "start": 5}}

objectives = {9: [{"measure": ["y", "alpha"],"goal":0,"weight":1},
                  {"measure": ["x", "alpha"],"goal":0,"weight":1},
                  {"measure": ["y", "std"],"goal":1,"weight":1},
                  {"measure": ["x", "std"],"goal":1,"weight":1}]}

matrixVariables = ebeam.gen_6d_gaussian(0,[1,.2,1,0.2,1,1],1000)
beam_dist = matrixVariables
test = beamOptimizer(line, vals, "SLSQP", matrixVariables, objectives, startPoint= starting)






# test = beamOptimizer(line, (1,3,5),1,1, [A,B], "Nelder-Mead", matrixVariables)


result = test.calc(plotProgress = True, beamlinePlotParams={"interval": 0.05})


# print("speed " + str(test.testSpeed(10)))
# evals, evalType = test.testFuncEval(10)
# print("evals " + str(evals) + ", " + evalType)
# print("iterations " + str(test.testFuncIt(10)))

line[1].current = result.x[0]
line[3].current = result.x[1]
line[5].current = result.x[1]
line[7].current = result.x[2]
schem.plotBeamPositionTransform(beam_dist, line, 0.05)

print("Var" + str(test.variablesToOptimize))
print("Current" + str(result.x))
print("MSE:" + str(result.fun))
