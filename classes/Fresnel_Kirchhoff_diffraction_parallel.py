"""
Fresnel-Kirchhoff diffraction
all distances in mm

"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)

import numpy as np
Pi=np.pi
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count


class FK_diffraction():
    def __init__(self,Ein,Xin,Yin,Xout,Yout,L,p,lam=1030*10**-6):
        self.Ein=Ein
        self.Xin=Xin
        self.Yin=Yin
        self.Xout=Xout
        self.Yout=Yout
        self.L=L
        self.p=p
        self.lam=lam
        self.k=2*Pi/self.lam
    
    Ein=None #electric field in the input plane
    Xin=None #starting plane coordinates
    Yin=None
    Xout=None #output plane coordinates
    Yout=None
    L=0 #Zout: distance bitween input and output planes
    Eout=None #electric field in the output plane
    lam=1030*10**-6 #wavelength
    p=None #multiprocess pool object
    lam=1030*10**-6 #wavelength
    k=None
    
    def propagate(self):
        Xin=self.Xin #assumed to be same as Yin
        # Yin=self.Yin
        lam=self.lam
        L=self.L
        dx=Xin[1]-Xin[0]
        Nx=len(Xin)
        XinM=Xin[:,None]*np.ones(Nx)
        YinM=np.ones(Nx)[:,None]*Xin
        p=self.p
        Ein=self.Ein
        Yout=self.Yout
        
        self.Eout=np.array(p.starmap(prop_x,[[Ein,XinM,YinM,x,Yout,L,lam,dx] for x in self.Xout]))
    
        
    
def prop_x(Ein,Xin,Yin,xout,Yout,L,lam,dx):
    """for parallelzation
    propagates 1 raw of data at fixed x position"""
    dy=dx
    Out=np.array([-1j/lam*dx*dy*np.sum(Ein*R_part(Xin,Yin,xout,y,L,2*Pi/lam)*
                                              Cos_part(Xin,Yin,xout,y,L))
                  for y in Yout])
    return Out
    
def R_part(xin,yin,xout,yout,z,k):
    """exp(i*k*R)/R"""
    R=((xin-xout)**2+(yin-yout)**2+z**2)**0.5
    # print(R)
    return np.exp(1j*k*R)/R

def Cos_part(xin,yin,xout,yout,z):
    """(1+cos(Rz))/2"""
    # Rvec=np.array([xin-xout,yin-yout,z])
    # Rvec=Rvec/np.sum(Rvec**2)**0.5 #normalization
    # Zv=np.array([0,0,1])
    # Cos=np.dot(Rvec,Zv)
    
    Cos=z/((xin-xout)**2+(yin-yout)**2+z**2)**0.5
    
    # return (1+Cos)/2
    return Cos
    
def F_2D(Ein,X,Kx,k,f,p):
    """fourier transform, f - focal length,
        p: Pool for parallel computations"""
    Nx=len(X)
    XinM=X[:,None]*np.ones(Nx)
    YinM=np.ones(Nx)[:,None]*X    
    return np.array(p.starmap(F_x,[[Ein,XinM,YinM,kx,Kx,k,f] for kx in Kx]))

def F_x(Ein,XinM,YinM,kx,Kx,k,f):
    """for parallelzation
    calculates 1 raw of data at fixed x position"""
    Out=np.array([np.sum(Ein*np.exp(-1j*k*kx*XinM/f-1j*k*ky*YinM/f))
                  for ky in Kx])
    return Out

def Fourier_2d_step(Ein,X,kx,k,f):
        """one focal point for the above function"""
        Nx=len(X)
        XinM=X[:,None]*np.ones(Nx)
        YinM=np.ones(Nx)[:,None]*X
        phase=np.exp(-1j*k*kx*XinM/f-1j*k*kx*YinM/f)
        return np.sum(Ein*phase)
    