"""
beam propagation
for hole reimaging (for IR-XUV separation)
Essentially consists of: 2D Fourier transform,-> applying an aperture and-> inverse 2D Fourier transform

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

#=========input parameters===========

W=3.5 #e**-2 radius of the beam. Be careful with increasing the beam size 
        #because the computation time might become too long
W_hole=W*0.2 #hole diameter (0 means no hole)
lam=1030*10**-6 #wavelength in mm
f=750 #in mm focal length of the spherical mirror
tilt=0 #in degree angle of incidence on the spherical mirror (=0.5 angle between incident and reflected beams)

Ra=300*10**-3 #radius of the aperture in the focus in mm
Nxin=151 #grid size (number of points in each dimention)

# print('e^-2 beam diameter: ', 2*W)
print('part of energy lost in the hole of the first mirror: ', 1-np.exp(-2*W_hole**2/W**2))

#========simulation parameters (automatically calculated based on the input parameters)=====

# W0=3.5/2**0.5
# Scaling=W/W0
DXin=W*3/3*2
# print('grid size: ',Nxin)
Xin=np.linspace(-DXin,DXin,Nxin)

#initiate the class
Fl=field()
Fl.GaussianBeam(Xin,W,lam=lam,W_hole=W_hole)
Fl.show_I('input beam')
# print(Fl.peak_fluence())
# print(Fl.diameter('e**-2'))

#if cycle is a wudu magic for parallel computations to work on windows
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
    if not tilt == 0: 
        Fl.add_tilt_mirror_aberration(tilt,2*f)
    Fl.Fourier(f,Ra*2.2,Nxin,p)
    Fl.show_I('focus')
    print('focus diameter (x, y):')
    print('e^-2 ',Fl.diameter('e**-2'))
    print('4sigma ',Fl.diameter('4sigma'))
    print('FWHM ',Fl.diameter('FWHM'))
    #apply aperture
    Fl.aperture(Ra)
    #back to near field
    Fl.Inv_Fourier(f,p)
    Fl.show_I('output beam')
    
    Portion=0.7 #how smaller is the hole of the second hole mirror
    print("part of energy in the full hole diameter: ", Fl.En_in_circle(W_hole))
    print('part of energy in the '+str(Portion)+' hole diameter: ',Fl.En_in_circle(Portion*W_hole))
    
    p.close()
    p.join()