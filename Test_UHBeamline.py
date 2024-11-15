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


pd.set_option('display.max_rows', None)
# Create beamline from Excel file
path2 = r"C:\Users\NielsB\cernbox\Hawaii University\Beam dynamics\FEL_sim"
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
line_E = beamtype.changeBeamType(beamlineUH, "electron", 45)



# ebeam
ebeam = beam()
beam_dist = ebeam.gen_6d_gaussian(0,[1,0.9,1,0.1,1,1],int(1e3))


schem.plotBeamPositionTransform(beam_dist, line_E, 0.1)