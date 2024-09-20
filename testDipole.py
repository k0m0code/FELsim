from ebeam import beam
from schematic import draw_beamline
from beamline import *
import numpy as np

ebeam = beam()
schem = draw_beamline()

I = 2.28
sec1 = driftLattice(0.05)
sec2 = dipole(length=0.0889, angle=1.5)
sec3 = driftLattice(0.050)
line = [sec1,sec2,sec3]

beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,1,1],10000)

schem.plotBeamPositionTransform(beam_dist, line, 0.001)
