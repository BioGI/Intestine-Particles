import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_file
#rc_file('rcFile')
from glob import glob

Vf_expectedData	= np.loadtxt('Plot_file_Vexpected.dat')
C_tot_Over_Cs	= Vf_expectedData[0,0]
Vp_tot		= Vf_expectedData[1,0] 
Np_tot          = int(Vf_expectedData[2,0])
dmin		= int(Vf_expectedData[3,0])
dmean		= int(Vf_expectedData[4,0])
dmax		= int(Vf_expectedData[5,0])
sigma		= int(Vf_expectedData[6,0])
D 	   	= Vf_expectedData[7:3000,0]
Vf_expected	= Vf_expectedData[7:3000,1]


nProcs = len(glob('Plot_file.dat'))
scalarData = []
for i in sorted(glob('Plot_file.dat')):
    scalarData.append(np.loadtxt(i))
Diameter    = scalarData[0][:100,0]
for i in range(nProcs):
    Vf_achieved = scalarData[i][:100,2] /Vp_tot
    Np          = scalarData[i][:100,3]

