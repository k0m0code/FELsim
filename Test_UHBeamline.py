#   Authors:
from pathlib import Path
from ebeam import beam
from schematic import draw_beamline
from excelElements import ExcelElements
from beamline import *
import sys
import time
import pandas as pd
from beamOptimizer import *
from AlgebraicOptimization import AlgebraicOpti
import sympy.plotting as plot
import sympy as sp

'''
Beam Properties initialization
Energy spread between 0.1 to 0.5 %
Bunch length 1 ps (1 deg) to 2 ps (2 deg) 
'''
ebeam = beam()
f = 2856 * (10 ** 6)
nb_particles = int(1e3)

bunch_spread = 1  # std in pico-second
tof_std = bunch_spread * (10 ** -9) * f  # (10 ** -3) (dToF / T)
energy_std = 1  # (10 ** -3) (dW / W)

x_std = 1  # (mm)
x_prime_std = 0.1  # (mrad)
y_std = 1  # (mm)
y_prime_std = 0.1  # (mrad)

beam_dist = ebeam.gen_6d_gaussian(0,
                                  [x_std,x_prime_std,y_std,y_prime_std,tof_std,energy_std],
                                  nb_particles)

'''
Import UH beamline lattice from an Excel file
Generate beamline() elements
'''

pd.set_option('display.max_rows', None)
# Create beamline from Excel file
path2 = r"C:\Users\NielsB\cernbox\Hawaii University\Beam dynamics\FELsim"
path1 = r"C:\Users\User\Documents\FELsim"
directory = Path(path2)
# file_path = directory / 'Beamline_elements.xlsx'
file_path = directory / 'Beamline_elements.xlsx'
excel = ExcelElements(file_path)
df = excel.get_dataframe()
#beamline
beamlineUH = excel.create_beamline()
# if len(beamline) >= 5:
#     beamline = beamline[:-5]
schem = draw_beamline()
beamtype = beamline()
line_UH = beamtype.changeBeamType(beamlineUH, "electron", 40)

print('Beamline nb of Elements: ' + str(len(line_UH)))
'''
Use optimization tools on segments of the UH beamline to the lattice values
'''
n = 66
truncatedA_line_optimization = line_UH[:-66]
truncatedB_line_optimization = line_UH[:-92]

j=0
for element in truncatedB_line_optimization:
    print(line_UH[j])
    j += 1

# Custom line
I = 3.56
sec_qpf = qpfLattice(current = I)
sec_qpd = qpdLattice(current = I)
sec_drift = driftLattice(0.50)
sec_dipole = dipole(length=0.0889, angle=1.5)
sec_wedge = dipole_wedge(length=0.01, angle=0, dipole_length=0.0889,dipole_angle=1.5)

line_test = [truncatedB_line_optimization,
            dipole_wedge(length=0.01, angle=0, dipole_length=0.0889,dipole_angle=1.5),
            dipole(length=0.0889, angle=1.5),
            dipole_wedge(length=0.01, angle=0, dipole_length=0.0889,dipole_angle=1.5),
            driftLattice(0.50),
            qpfLattice(current = I),
            driftLattice(0.50),
            dipole_wedge(length=0.01, angle=0, dipole_length=0.0889, dipole_angle=1.5),
            dipole(length=0.0889, angle=1.5),
            dipole_wedge(length=0.01, angle=0, dipole_length=0.0889, dipole_angle=1.5),
            driftLattice(0.50),
            qpdLattice(current = I),
            driftLattice(0.50),
            qpfLattice(current = I),
            driftLattice(0.50),
            qpdLattice(current = I),
            driftLattice(0.50)]

z_dpw1 =1.697921
'''
Plotting results
'''

twiss_aggregated_df = schem.plotBeamPositionTransform(beam_dist, truncatedB_line_optimization, 1.0)

print(twiss_aggregated_df)

#schem.plotBeamPositionTransform(beam_dist, line_test, 0.1)
