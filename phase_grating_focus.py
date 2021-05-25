"""
Created on Fri Mar 19 15:27:30 2021

@author: Slawa
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)
sys.path.append(Path+'//classes')
import numpy as np
Pi=np.pi
import matplotlib.pyplot as plt
from field_class import field
from multiprocessing import Pool, cpu_count

#=====phase grating parameters====
period=2 #grating period in mm
shift=period*0 #grating 


#=========input parameters===========

W=5 #e**-2 radius of the beam.
lam=800*10**-6 #wavelength in mm
f=700 # focal length in mm


Ra=200*10**-3 #radius of the image in focus in mm

Nxin=151 #grid size (number of points in each dimension)
# Nxin=41

#========simulation parameters (automatically calculated based on the input parameters)=====

DXin=W*2
Xin=np.linspace(-DXin,DXin,Nxin)

#initiate the class
Fl=field()
Fl.GaussianBeam(Xin,W,lam=lam)
Fl.add_phase_grating(period,shift)
Fl.show_I('input beam, intensity')
Fl.show_Phase('input beam, phase')

if __name__ == '__main__':
    #initiate parallel computing
    Ncpu=cpu_count()
    if Ncpu > 1:
        Npr=Ncpu-1 #number of parallel processes in the Pool
    else:
        Npr=1
    p=Pool(Npr)
    
    #==========computation============
    
    #far field
    Fl.Fourier(f,Ra*2,Nxin,p) #propagate
    Fl.show_I('focus') #show the plot
    
    p.close()
    p.join()

# L=f*3 #distance from the phase plate
# W0=3.5/2**0.5
# Scaling=W/W0
# DXin=W*3/3*2
# Nxin=int(46*Scaling**2*3000/L)+1
# print('grid size: ',Nxin)
# Xin=np.linspace(-DXin,DXin,Nxin)

# #initiate the class
# Fl=field()
# Fl.GaussianBeam(Xin,W,lam=lam)
# Fl.add_phase_grating(period,shift)
# Fl.show_I('input beam, intensity')
# Fl.show_Phase('input beam, phase')

# if __name__ == '__main__':
#     #initiate parallel computing
#     Ncpu=cpu_count()
#     if Ncpu > 1:
#         Npr=Ncpu-1 #number of parallel processes in the Pool
#     else:
#         Npr=1
#     p=Pool(Npr)
    
#     Fl.FK_propagation(L,p)
#     Fl.show_I('after phase grating on: ',L,' mm distance')