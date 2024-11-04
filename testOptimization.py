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

vals = {2: ["L","length", lambda num:num],
        1: ["I", "current", lambda num:num],
        3: ["I", "current", lambda num:num],
        5: ["I", "current", lambda num:num*2],
        7: ["I", "current", lambda num:num]}

starting = {"I": {"bounds": (0,10), "start": 5},
            "L": {"bounds": (0,1), "start": 0.1}}

objectives = {9:[{"measure": ["y", "std"],"goal":1,"weight":1},
                 {"measure": ["x", "std"],"goal":1,"weight":1}]}

matrixVariables = ebeam.gen_6d_gaussian(0,[1,.2,1,0.2,1,1],1000)
beam_dist = matrixVariables
test = beamOptimizer(pBeam, matrixVariables)

result = test.calc("Nelder-Mead", vals, starting, objectives, plotProgress = True, plotBeam= True, printResults=True)







# test = beamOptimizer(line, (1,3,5),1,1, [A,B], "Nelder-Mead", matrixVariables)


