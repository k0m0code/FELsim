#   Authors:
from pathlib import Path 
from beamline import lattice
from ebeam import beam
from beamline import *
from schematic import *
import random
import numpy as np

# Define the directory and file path using pathlib
directory = Path("C:/Users/User/Documents/FELsim")
file_path = directory / 'Beamline_elements.xlsx'

# Load beamline lattice and electron beam properties
beamline = lattice()
ebeam = beam()


# beamline.load_excel_lattice(file_path)
# print(beamline.nomenclatures)
# test = beamline.find_element_by_position(4.19)
# print(beamline.positions)
# print(test)

#print(beamline)
#print(beamline.positions.head())  # Print the first few rows of the positions DataFrame


# z_start_array = beamline.positions['z_sta'].values
# element = beamline.nomenclatures

#dist = ebeam.gen_6d_gaussian(5,1.6)
#lista = list[:,0]
#listb = list[:,1]
#listc = list[:,2]
#listd = list[:,3]
# liste = list[:,4]
# listf = list[:,5]

# test = [["drift", "QPF", "drift", "QPF"],
#          [500,88.9,900,88.9]]
# # ebeam.plotBeamPositionTransform(list, test,150)

sec1 = driftLattice(200)
sec2 = qpfLattice(current = 32)
sec3 = driftLattice(100)
sec4 = qpdLattice(current = 32)
sec5 = driftLattice(100)
sec6 = qpfLattice(current = 32)
sec7 = driftLattice(100)
sec8 = qpdLattice(current = 32)
sec9 = driftLattice(200)
beamline_lattice = [sec1,sec2,sec3,sec4,sec5,sec6,sec7,sec8,sec9]
beam_dist = ebeam.gen_6d_gaussian(0,[1,.1,1,0.1,1,1],1000)



# ebeam.particles_in_ellipse(beam_dist)
# ebeam.cal_twiss(beam_dist)
# ebeam.plot_6d(beam_dist, 'bruh')

schem = draw_beamline()
# tup = ebeam.getXYZ(beam_dist)

# schem.driftTransformScatter(beam_dist,70,True)
# schem.plotBeamPositionTransform(beam_dist, beamline_lattice, 5, saveData=True)
# schem.plotBeamPositionTransform(beam_dist, beamline_lattice, 500, plot_z = (400, 0))
shape1 = {"shape": "circle", "radius": 5, "origin": (0,5)}
shape = {"shape": "bruh", "length": 5, "width": 10, "origin": (10,-4)}
schem.plotBeamPositionTransform(beam_dist, beamline_lattice, 10, True, shape = shape)

# print(beam_dist)
# bru = sec1.useMatrice(beam_dist,100)
# print(bru)

# test = [["drift","drift","drift"],[3,3,1.6]]
# ebeam2.plotBeamPositionTransform(list, test, 2)

# test = np.array([[1,2,3,4,5,6],[1,1,1,1,1,1]])
# result = ebeam.getDriftMatrice(test, length = 5)
# print(result)


# ebeam.plot_6d(lista, listb, listc,listd, 
#               liste, listf)
# newarray = []
# for i in list:
#     newarray.append(ebeam.getDriftMatrice(i))
# transformedMatrix = np.array(newarray)
# ebeam.plot_6d(transformedMatrix[:,0],transformedMatrix[:,1],transformedMatrix[:,2],transformedMatrix[:,3],transformedMatrix[:,4],transformedMatrix[:,5])
# ebeam.plotDriftTransform(lista,listc,listb,listd)
