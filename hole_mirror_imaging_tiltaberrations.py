"""
beam propagation
Fresnel-Kirchhoff integral for hole reimaging (for IR-XUV separation)

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
tilt=3 #in degree angle of incidence on the spherical mirror (=0.5 angle between incident and reflected beams)
L1=1500 #distance between the first hole mirror and the spherical mirror in mm
L2=1/(1/f-1/L1) #distance to the second hole mirror from the imaging condition

# print('e^-2 diameter: ', 2*W)
print('part of energy lost in the hole of the first mirror: ', 1-np.exp(-2*W_hole**2/W**2))
print('L2: ',L2)

#========simulation parameters (automatically calculated based on the input parameters)=====

W0=3.5/2**0.5
Scaling=W/W0
DXin=W*3/3*2
Nxin=int(46*Scaling**2*3000/min([L1,L2]))+1
print('grid size: ',Nxin)
Xin=np.linspace(-DXin,DXin,Nxin)

#initiate the class
Fl=field()
Fl.GaussianBeam(Xin,W,lam=lam,W_hole=W_hole)
Fl.show_I('input beam')
# print(Fl.peak_fluence()) #peak fluence in the beam
# print(Fl.diameter('e**-2')) #beam size

#==========simulation=============
#if cycle is a wudu magic for parallel computations to work on windows
if __name__ == '__main__':
    #initiate parallel computing
    Ncpu=cpu_count()
    if Ncpu > 1:
        Npr=Ncpu-1 #number of parallel processes in the Pool
    else:
        Npr=1
    p=Pool(Npr)
    
    #=============Fresnel-Kirchhoff propagation============
    #propagate the the spherical mirror
    Fl.aperture(DXin)
    Fl.FK_propagation(L1,p)
    Fl.show_I('on spherical mirror')
    
    #add spherical mirror
    Fl.add_lens(f)
    #add tilt aberrations (remove if no aberrations are needed)
    if not tilt == 0:
        Fl.add_tilt_mirror_aberration(tilt,2*f)
    
    #propagate to the second hole mirror
    Fl.FK_propagation(L2,p)
    Fl.show_I('output beam')
    
    #calculate the energy content leaking through the hoe of the second hole mirror
    Portion=0.7 #how smaller is the hole of the second hole mirror
    print("part of energy in the full hole1 diameter: ", Fl.En_in_circle(0.98*W_hole))
    print('part of energy in the '+str(Portion)+' hole1 diameter: ',Fl.En_in_circle(Portion*W_hole))
