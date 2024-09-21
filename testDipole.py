from ebeam import beam
from schematic import draw_beamline
from beamline import *
import numpy as np

ebeam = beam()
schem = draw_beamline()

I = 2.28
sec1 = driftLattice(0.01)
sec2 = dipole_wedge(length=0.01, angle=15, dipole_length=0.0889,dipole_angle=15)
sec3 = dipole(length=0.0889, angle=15)
#sec4 = dipole_wedge(length=0.01, angle=15, dipole_length=0.0889,dipole_angle=15)
sec12 = driftLattice(0.01)
line = [sec1,sec2,sec12]

beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,1,1],10000)

schem.plotBeamPositionTransform(beam_dist, line, 0.01)
