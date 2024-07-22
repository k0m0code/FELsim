from ebeam import beam
import random
import numpy as np

ebeam = beam()


# '''
# Testing matrix multiplcation for driftspace
# '''
# E = 35  # Kinetic energy (MeV/c^2)
# E0 = 0.51099
# length = 3
# gamma = (1 + (E/E0))
# list1 = np.array([-2.61479201,
#                  2.86977149,
#                 -5.7357132,
#                  2.94133543,
#                  0.86818946,
#                 -0.47795909])
# list1 = list1.T
# print(list1)
# print(driftMatrice)
# print(np.matmul(driftMatrice,list1))

list = ebeam.gen_6d_gaussian(0,random.randint(1,30), 100)
beamline = [["drift","drift",'drift','drift'],[1,1,3,5]]
ebeam.plotBeamPositionTransform(list,beamline)