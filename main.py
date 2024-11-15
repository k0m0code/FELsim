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
beam_dist = ebeam.gen_6d_gaussian(0,[1,0.1,1,0.1,1,1],10000)


schem.plotBeamPositionTransform(beam_dist, line_E, 10)

'''
create beamline, get sigma i
'''
alg = AlgebraicOpti()
sig = alg.getDistSigmai(beam_dist)

'''
Initial conditions, sigma_i for horizontal and vertical planes:
'''
print("Initial h-plane sigma_i[0] (eps * beta):"+str(sig[0]))
print("Initial h-plane sigma_i[1] (-eps * alpha):"+str(sig[1]))
print("Initial h-plane sigma_i[6] (-eps * alpha):"+str(sig[6]))
print("Initial h-plane sigma_i[7] (eps * gamma):"+str(sig[7]))
print("Initial v-plane sigma_i[12] (eps * beta):"+str(sig[14]))
print("Initial v-plane sigma_i[13] (-eps * alpha):"+str(sig[15]))
print("Initial v-plane sigma_i[18] (-eps * alpha):"+str(sig[20]))
print("Initial v-plane sigma_i[19] (eps * gamma):"+str(sig[21]))
'''
We would like to compare this values with sigma_f, to try to have sigma_f = sigma_i for the h- and v-plane
'''
print(sig)

I = 3.56  # result obtained from testOptimization.py...
sec1 = driftLattice(0.5)
sec2 = qpfLattice(current = I)
sec3 = driftLattice(0.25)
sec4 = qpdLattice(current = I)
sec5 = driftLattice(0.25)
sec6 = qpfLattice(current = I)
sec7 = driftLattice(0.25)
sec8 = qpdLattice(current = I)
sec9 = driftLattice(0.50)
line = [sec1,sec2,sec3,sec4,sec5,sec6,sec7,sec8,sec9]

beamtype = beamline()
line_E = beamtype.changeBeamType(line, "electron", 55)

'''
plotting
'''
# schem.plotBeamPositionTransform(beam_dist,line_E,0.1, spacing = 5)


'''
create x values to optimize
{segment parameter: variable name}
'''
xvals = {
         1: {"current": "I"},
         3: {"current": "I"},
         5: {"current": "I"},
         7: {"current": "I"}
        }



mAr = alg.getM(line, xvals)
'''
def objectives
'''
#epsilon, alpha, beta, gamma
yObj = {'x': [0,0,0,0],'y': [0.999,0,0.02,0.23],'z': [0.9, 0.003, 1, 0.2]}
finm = alg.findObj(line_E, xvals, yObj)
print(finm[15])
plot.plot(finm[15])
I = sp.symbols("I", real = True)
print(sp.solveset(finm[15], I, domain=sp.S.Reals))

'''
print out the M matrice of the beamline
'''
# for row in mAr.tolist():
#     print(row)

'''
print values for sigmaF, and plot each function
'''
sigmaf = alg.getSigmaF(mAr,sig)

# for row in sigmaf.tolist():
#     print("sigmaf:"+str(row) + "\n")

nb_points = 1000  # number of point used in the simpy


'''
This figures are for testing purposes.
The figure that is interesting should be the objective function as a function of the variable,
for instance here, F(I) = abs(sigma_i(1,1) - sigma_f(1,1)) + abs(sigma_i(5,5)-sigma_f(5,5)), 
'''
p1 = plot.plot(sigmaf[0], nb_of_points=nb_points, line_color='red', show=False) # find I so that this equals ~ 1
p2 = plot.plot(sigmaf[1], nb_of_points=nb_points,line_color='blue', show=False) # find I so that this equals ~ 0
p3 = plot.plot(sigmaf[7], nb_of_points=nb_points,line_color='green', show=False) # find I so that this equals ~ 0
sig_i_0 = plot.plot(sig[0], line_color='red', show=False)
sig_i_1 = plot.plot(sig[1], line_color='blue', show=False)
sig_i_7 = plot.plot(sig[7], line_color='green', show=False)

p1.append(p2[0])
p1.append(p3[0])
p1.append(sig_i_0[0])
p1.append(sig_i_1[0])
p1.append(sig_i_7[0])

p1.show()




'''
optimizization example
'''
# optimize = beamOptimizer(beamline, beam_dist)
# segmentVar = {1: ["I", "current", lambda num:num],
#               3: ["B", "current", lambda num:num]}
# objectives = {4: [{"measure": ["y", "alpha"],"goal":0,"weight":1},
#                   {"measure": ["x", "std"],"goal":1,"weight":1}]}
# starting = {"I": {"bounds": (0,10), "start" :5}}
# optimize.calc("COBYQA", segmentVar, starting, objectives, True, True, True)



'''
Shape example
'''
# shape1 = {"shape": "circle", "radius": 5, "origin": (0,5)}
# shape = {"shape": "rectangle", "length": 200, "width": 500, "origin": (10,-4)}
# schem.plotBeamPositionTransform(beam_dist, beamline, 10, shape = shape)


'''
test
'''
# x = smp.symbols("x")
# eq = x**2
# xp = plot.plot(eq, show = False)
# xp.show()
# eq2 = eq - 5
# xp2 = plot.plot(eq2, show = False)
# xp.append(xp2[0])
# xp.show()
