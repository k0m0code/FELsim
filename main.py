#   Authors:
from pathlib import Path
from ebeam import beam
from schematic import draw_beamline
from excelElements import ExcelElements
import numpy as np

# # Create beamline from Excel file
# path2 = r"C:\Users\NielsB\cernbox\Hawaii University\Beam dynamics\FEL_sim"
# path1 = r"C:\Users\User\Documents\FELsim"
# directory = Path(path1)
# file_path = directory / 'Beamline_elements.xlsx'
# excel = ExcelElements(file_path)
# df = excel.get_dataframe()
# print(df)


#beamline
# beamline = excel.create_beamline()
# if len(beamline) >= 5:
#     beamline = beamline[:-5]
schem = draw_beamline()

print(10**(np.ceil(np.log10(0.000012342))))


# ebeam
ebeam = beam()
beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,1,1],1000)
tup = (ebeam.cal_twiss, 1)



dist_avg, dist_cov, twiss =tup[0](beam_dist, ddof=1)
print(twiss)
print(twiss.loc['x'])
print(twiss.loc["x"].loc[r"$\alpha$"])

# schem.plotBeamPositionTransform(beam_dist, beamline, 0.2)



# '''
# Shape example
# '''

# shape1 = {"shape": "circle", "radius": 5, "origin": (0,5)}
# shape = {"shape": "rectangle", "length": 200, "width": 500, "origin": (10,-4)}
# schem.plotBeamPositionTransform(beam_dist, beamline, 10, shape = shape)
