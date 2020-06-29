"""
Fresnel-Kirchhoff diffraction
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)

import numpy as np
Pi=np.pi
import matplotlib.pyplot as plt
from Fresnel_Kirchhoff_diffraction_parallel import FK_diffraction
from multiprocessing import Pool, cpu_count


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
        
def lens(x,y,f,lam=1030*10**-6):
    """phase of a thin lens with focus f"""
    k=2*Pi/lam
    return np.exp(-1j*k/2/f*(x**2+y**2))

def x_out(xin,L,lam=1030*10**-6):
    k=2*Pi/lam
    dxin=xin[1]-xin[0]
    DXout=0.5/dxin/k*L*2*Pi
    return np.linspace(-DXout,DXout,len(xin))

DXin=W*3/3*2
Nxin=int(46*Scaling**2*3000/L)+1
Xin=np.linspace(-DXin,DXin,Nxin)
dxin=Xin[1]-Xin[0]

print(2*DXin,' max shift', 2*Pi/dxin/k*L, ' shift step', 2*Pi/2/DXin/k*L)

Ein=np.array([[E_Gaus(x,y,W,R,lam,W_hole) for y in Xin] for x in Xin])

plt.pcolormesh(Xin,Xin,np.abs(Ein)**2)
plt.show()

# DXout=W*3
DXout=0.5/dxin/k*L*2*Pi
# Nxout=Nxin
# Nxout=30+1
Xout=x_out(Xin,L,lam)

print('size in ',2*DXin,'size out', 2*DXout,' max shift', 2*Pi/dxin/k*L, ' shift step', 2*Pi/2/DXin/k*L)

#if cycle is a wudu magic for parallel computations to work in windows
if __name__ == '__main__':
    #initiate parallel computing
    Ncpu=cpu_count()
    if Ncpu > 1:
        Npr=Ncpu-1 #number of parallel processes in the Pool
    else:
        Npr=1
    p=Pool(Npr)
    
    #propagarte to lens
    Dif=FK_diffraction(Ein,Xin,Xin,Xout,Xout,L,p)
    Dif.propagate()
    
    plt.pcolormesh(Xout,Xout,np.abs(Dif.Eout)**2)
    plt.show()
    
    #add lens
    f=750
    Ein2=Dif.Eout*np.array([[lens(x,y,f,lam) for y in Xout] for x in Xout])
    
    L2=1500
    Xin2=Xout
    Xout2=x_out(Xin2,L2,lam)
    
    Dif2=FK_diffraction(Ein2,Xin2,Xin2,Xout2,Xout2,L2,p)
    Dif2.propagate()
    
    plt.pcolormesh(Xout2,Xout2,np.abs(Dif2.Eout)**2)
    plt.show()
    
    
    
    