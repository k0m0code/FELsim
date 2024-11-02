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
from AlgebriacOptimization import AlgebriacOpti

pd.set_option('display.max_rows', None)
# Create beamline from Excel file
# path2 = r"C:\Users\NielsB\cernbox\Hawaii University\Beam dynamics\FEL_sim"
path1 = r"C:\Users\User\Documents\FELsim"
directory = Path(path1)
file_path = directory / 'Beamline_elements(1).xlsx'
excel = ExcelElements(file_path)
df = excel.get_dataframe()


#beamline
beamline = excel.create_beamline()
# if len(beamline) >= 5:
#     beamline = beamline[:-5]
schem = draw_beamline()


seg = driftLattice(1)
print(seg.getSymbolicMatrice(**{"length": "L"}))


# ebeam
ebeam = beam()
beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,1,1],1000)

alg = AlgebriacOpti()
sig = alg.getSigmai(beam_dist)


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
line = [sec1,sec2]

xvals = {0:{"length":"BRUHH"}}
ar = alg.getM(line, xvals)
for row in ar.tolist():
    print(row)


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
