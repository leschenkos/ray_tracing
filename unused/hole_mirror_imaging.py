"""
beam propagation
main part is hole reimaging (for IR-XUV separation)

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


Scaling=2**0.5
W=1.987*Scaling*1.5 #radius e-2
W=3.5 #e**-2 radius of the beam. Be careful with increasing the beam size 
        #because the computation time might become too long
W_hole=W*0.25*1 #hole diameter (0 means no hole)
print('e-2 radius', W)
print('part of energy lost in the hole of the first mirror: ', 1-np.exp(-2*W_hole**2/W**2))
lam=1030*10**-6 #wavelength in mm
L=1500 #distance in mm
DXin=W*3/3*2
Nxin=int(46*Scaling**2*3000/L)+1
# print(Nxin)
Xin=np.linspace(-DXin,DXin,Nxin)

Fl=field()
Fl.GaussianBeam(Xin,W,lam=lam,W_hole=W_hole)
Fl.show_I()
print(Fl.peak_fluence())
# print(Fl.diameter('e**-2'))

#if cycle is a wudu magic for parallel computations to work in windows
if __name__ == '__main__':
    #initiate parallel computing
    Ncpu=cpu_count()
    if Ncpu > 1:
        Npr=Ncpu-1 #number of parallel processes in the Pool
    else:
        Npr=1
    p=Pool(Npr)
    
    #=============FK propagation and mirror tilt aberrations
    # print(p._processes)
    # Fl.aperture(DXin)
    # Fl.FK_propagation(L,p)
    # Fl.show_I('on spherical mirror')
    # print(Fl.En_in_circle(0.9*W_hole))
    # print(Fl.peak_fluence())
    
    # tilt=10
    # f=750
    
    # Fl.add_lens(f)
    # Fl.add_tilt_mirror_aberration(tilt,2*f)
    
    # Fl.FK_propagation(L,p)
    # Fl.show_I()
    # print(Fl.En_in_circle(0.98*W_hole))
    
    #==========far field============
    f=750 #focal length of the spherical mirror in mm
    Ra=325*10**-3 #radius of the aperture in the focus in mm
    tilt=4 #angle of incidence on the spherical mirror in degree (=0.5 angle between incident and reflected beams)
    Fl.add_tilt_mirror_aberration(tilt,2*f)
    Fl.Fourier(f,Ra*2.2,Nxin,p)
    Fl.show_I('focus')
    print('focus diameter (x, y):')
    print('e^-2 ',Fl.diameter('e**-2'))
    print('4sigma ',Fl.diameter('4sigma'))
    print('FWHM ',Fl.diameter('FWHM'))
    Fl.aperture(Ra)
    Fl.Inv_Fourier(f,p)
    Fl.show_I('image on the output')
    
    Portion=0.7 #how smaller is the hole of the second hole mirror
    print("part of energy in the full hole diameter: ", Fl.En_in_circle(W_hole))
    print('part of energy in the '+str(Portion)+' hole diameter: ',Fl.En_in_circle(Portion*W_hole))
    
    p.close()
    p.join()