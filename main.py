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

pd.set_option('display.max_rows', None)
# Create beamline from Excel file
# path2 = r"C:\Users\NielsB\cernbox\Hawaii University\Beam dynamics\FEL_sim"
path1 = r"C:\Users\User\Documents\FELsim"
directory = Path(path1)
file_path = directory / 'Beamline_elements(1).xlsx'
excel = ExcelElements(file_path)
df = excel.get_dataframe()
print(df)


#beamline
beamline = excel.create_beamline()
# if len(beamline) >= 5:
#     beamline = beamline[:-5]
schem = draw_beamline()

test = qpdLattice(2.4)
mat = test.getSymbolicMatrice()
equation = mat[0]
print(equation.subs({I:2.3,l:1}))

# ebeam
ebeam = beam()
beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,1,1],1000)

# schem.plotBeamPositionTransform(beam_dist,beamline,0.1, spacing = 5)

# optimize = beamOptimizer(beamline, beam_dist)
# segmentVar = {1: ["I", "current", lambda num:num],
#               3: ["B", "current", lambda num:num]}
# objectives = {4: [{"measure": ["y", "alpha"],"goal":0,"weight":1},
#                   {"measure": ["x", "std"],"goal":1,"weight":1}]}
# starting = {"I": {"bounds": (0,10), "start" :5}}
# optimize.calc("COBYQA", segmentVar, starting, objectives, True, True, True)



# '''
# Shape example
# '''

# shape1 = {"shape": "circle", "radius": 5, "origin": (0,5)}
# shape = {"shape": "rectangle", "length": 200, "width": 500, "origin": (10,-4)}
# schem.plotBeamPositionTransform(beam_dist, beamline, 10, shape = shape)
