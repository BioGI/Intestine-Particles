import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import commands
from matplotlib.patches import Ellipse, Polygon
import numpy as np
from matplotlib import rc_file
from pylab import loadtxt
import pickle, sys, os
matplotlib.rcParams.update({'font.size': 10})

title_font = {'fontname':'Arial', 'size':'18', 'color':'black', 'weight':'normal','verticalalignment':'bottom'} 
axis_font = {'fontname':'Arial', 'size':'14'}

os.chdir('../Case-0')
import Plot0 as Case0
#os.chdir('../Case-1')
#import Plot1 as Case1
#os.chdir('../Case-3')
#import Plot3 as Case3
#os.chdir('../Case-5')
#import Plot5 as Case5
#os.chdir('../Case-6')
#import Plot6 as Case6
#os.chdir('../Case-7')
#import Plot7 as Case7
os.chdir('../Plots')

#============================================================================================================
fig = plt.figure()
plt.axes([0.13,0.13,0.97-0.13,1.0-0.23])
plt.bar(Case0.Diameter-2.5, Case0.Np, 5,  label= "")
plt.grid(b=True, which='major',  axis='y', linewidth=1)
plt.title(r'$\frac{C_{tot}}{C_s}$=%s, $D_{min}$= %s$\mu m$,  $D_{max}$= %s$\mu m$,  $\sigma$= %s$\mu m$,  $N_P$=%s'  %( Case0.C_tot_Over_Cs, Case0.dmin, Case0.dmax, Case0.sigma, Case0.Np_tot), **title_font)
plt.xlim(0,200)
plt.xlabel(r'Diameter ($\mu m$)',**axis_font)
plt.ylabel('Number of Particles in each bin',**axis_font)
plt.legend(loc=0)
plt.savefig('Np.png')


#============================================================================================================
fig = plt.figure()
plt.axes([0.13,0.13,0.97-0.13,1.0-0.23])
plt.plot(Case0.D, Case0.Vf_expected, color='blue', label= "Expecetd") 
plt.bar(Case0.Diameter-2.5, Case0.Vf_achieved, 5, color='red', label= "Achieved")
plt.title(r'$\frac{C_{tot}}{C_s}$=%s, $D_{min}$= %s$\mu m$,  $D_{max}$= %s$\mu m$,  $\sigma$= %s$\mu m$,  $N_P$=%s'  %( Case0.C_tot_Over_Cs, Case0.dmin, Case0.dmax, Case0.sigma, Case0.Np_tot), **title_font)
plt.xlim(0,200)
plt.xlabel(r'Diameter ($\mu m$)',**axis_font)
plt.ylabel(r'PDF ($1/ \mu m)$ ',**axis_font)
plt.legend(loc=0)
plt.savefig('Vf.png')

#============================================================================================================
commands.getoutput("convert -border 1x2 +append Vf.png Np.png Particles.png")
commands.getoutput("chromium  Particles.png")

