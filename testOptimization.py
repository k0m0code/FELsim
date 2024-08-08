from beamOptimizer import beamOptimizer
from ebeam import beam
from schematic import draw_beamline
from beamline import *

ebeam = beam()
beam_dist = ebeam.gen_6d_gaussian(0,[1,.1,1,0.1,1,1],1000)
sec1 = driftLattice(200)
sec2 = qpfLattice(current = 32)
sec3 = driftLattice(100)
sec4 = qpdLattice(current = 32)
line = [sec1,sec2,sec3,sec4]
test = beamOptimizer(beam_dist,line, (0,4),25,50)

result = test.calc([1,1])
print(result.x)
print(result.fun)


schem = draw_beamline()

line[1].current = result.x[0]
line[2].current = result.x[1]
schem.plotBeamPositionTransform(beam_dist, line, 10)

