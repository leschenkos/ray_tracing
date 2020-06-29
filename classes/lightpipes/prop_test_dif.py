# -*- coding: utf-8 -*-
"""
Created on Thu May 28 17:00:33 2020

@author: Slawa
"""
import numpy as np
Pi=np.pi
import matplotlib.pyplot as plt
from LightPipes import *
from LightPipes.field import Field

#input parameters
Scaling=2**0.5
W=1.987*Scaling*1.5 #radius e-2
W_hole=W*0.1
print('e-2 radius', W)
lam=1030*10**-6 #wavelength
k=2*Pi/lam
R=10**15 #wavefront curvature
L=1500 #distance


def E_Gaus(x,y,W,R,lam=1030*10**-6,Whole=0):
    """field of a gaussian beam"""
    r=(x**2+y**2)**0.5
    k=2*Pi/lam
    if Whole==0:
        return np.exp(-r**2/W**2-1j*k*r**2/2/R)
    else:
        if r <= Whole:
            return 0
        else:
            return np.exp(-r**2/W**2-1j*k*r**2/2/R)

DXin=W*3/3*2
Nxin=int(46*Scaling**2*3000/L)+1
Xin=np.linspace(-DXin,DXin,Nxin)
dxin=Xin[1]-Xin[0]

print(2*DXin,' max shift', 2*Pi/dxin/k*L, ' shift step', 2*Pi/2/DXin/k*L)

Ein=np.array([[E_Gaus(x,y,W,R,lam,W_hole) for y in Xin] for x in Xin])

plt.pcolormesh(Xin,Xin,np.abs(Ein)**2)
plt.show()

wavelength=lam*10**3*um
size=DXin*2*mm
N=Nxin

F=Field(Ein,size,wavelength,N)
# F2=Forward(L*mm,size,N,F)
F2=Fresnel(L*mm,F)
# F2=Forvard(L*mm,F)
I=Intensity(2,F2)

plt.imshow(I,cmap='jet')
plt.show()