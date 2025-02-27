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
f = 2856 * (10 ** 6)  # Accelerator RF frequency (Hz)

nb_particles = int(1e4)

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
path3 = r"/Users/christiankomo/Desktop/Documents/FELsim"
path2 = r"C:\Users\NielsB\cernbox\Hawaii University\Beam dynamics\FELsim"
path1 = r"C:\Users\User\Documents\FELsim"
directory = Path(path3)
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

segments = 39
line = line_UH[:segments]
opti = beamOptimizer(line, beam_dist)
# schem.plotBeamPositionTransform(beam_dist, line,0.01)


# variables = {1: ["I", "current", lambda num:num],
#              3: ["I2", "current", lambda num:num]}
# startPoint = {"I": {"bounds": (0,10), "start": 1},
#            "I2": {"bounds": (0, 10), "start": 1}}
# objectives = {8: [{"measure": ["x", "alpha"], "goal": -5, "weight":1},]
#                #  {"measure": ["y", "beta"], "goal": 1, "weight":1}],
#                # 9:[{"measure": ["x", "alpha"], "goal": 0, "weight":1},
#                #  {"measure": ["y", "alpha"], "goal": 0, "weight":1}]
#              }
# result = opti.calc("Nelder-Mead", variables, startPoint, objectives, plotBeam= True, printResults=True)
line[1].current =  0.9989681933
line[3].current = 1.044851479

# variables = {10: ["I", "current", lambda num:num]}
# startPoint = {"I": {"bounds": (0,10), "start": 1}}

# objectives = {10: [{"measure": ["x", "alpha"], "goal": 5, "weight": 1},
#                    {"measure": ["y", "alpha"], "goal": -5, "weight": 20}]}
#
# result = opti.calc("Nelder-Mead", variables, startPoint, objectives, plotBeam= True, printResults=True)
line[10].current = 3.5

# variables = {16: ["I", "current", lambda num:num],
#              18: ["I2", "current", lambda num:num],
#              20: ["I", "current", lambda num:num]}
# startPoint = {"I": {"bounds": (0,10), "start": 2},
#               "I2": {"bounds": (0, 10), "start": 2}}

# objectives = {24: [{"measure": ["x", "envelope"], "goal": 0, "weight": 1},
#                    {"measure": ["y", "alpha"], "goal": 5, "weight": 10}],
#               16: [{"measure": ["y", "envelope"], "goal": 1.4, "weight": 1}],
#               20: [{"measure": ["y", "envelope"], "goal": 2, "weight": 1}]}
#
# result = opti.calc("Nelder-Mead", variables, startPoint, objectives, plotBeam= True, printResults=True, plotProgress=True)
line[16].current =  2.882858241
line[18].current =  5.108214683
line[20].current = 2.882858241

schem.plotBeamPositionTransform(beam_dist, line,0.01)

variables = {27: ["I", "current", lambda num:num]}
startPoint = {"I": {"bounds": (0,10), "start": 1}}

objectives = {
               # 32: [{"measure": ["y", "envelope"], "goal": 0.7, "weight": 1}],
               28: [{"measure": ["x", "envelope"], "goal": 0, "weight": 10}],
               27: [{"measure": ["y", "alpha"], "goal": -2.5, "weight": 1}]
              }

result = opti.calc("Nelder-Mead", variables, startPoint, objectives, plotBeam= True, printResults=True, plotProgress=True)

schem.plotBeamPositionTransform(beam_dist, line,0.01)
for i in line: print(i)



'''
ALGEBRAIC OPTIMIZATION WORKPLACE BELOW
'''
# testLine = line_UH[16:21]
# alg = AlgebraicOpti()
# xvar = {
#          0: {"current": "I2"},
#          2: {"current": "I"},
#          4: {"current": "I2"}
# }
# alg.BIVARIATE_SEARCH_RANGE = 5
# finm = alg.findSymmetricObjective(testLine, xvar, beam_dist, plotBeam= [2,3])
# print(finm[0,0])
