#   Authors:
from pathlib import Path
from ebeam import beam
from schematic import draw_beamline
from excelElements import ExcelElements
from beamline import lattice

# Create beamline from Excel file
path2 = r"C:\Users\NielsB\cernbox\Hawaii University\Beam dynamics\FEL_sim"
path1 = "C:/Users/User/Documents/FELsim"
directory = Path(path2)
file_path = directory / 'Beamline_elements.xlsx'
excel = ExcelElements(file_path)
df = excel.get_dataframe()
print(df)
#beamlattice = lattice.changeE(E=70)

beamline = excel.create_beamline()



#if len(beamline) >= 5:
#    beamline = beamline[:-5]
schem = draw_beamline()


# ebeam
ebeam = beam()
beam_dist = ebeam.gen_6d_gaussian(0,[1,1,1,1,0.1,450],1000)

schem.plotBeamPositionTransform(beam_dist, beamline, 0.4, defineLim=True, saveData=False, shape = {}, plot=True)