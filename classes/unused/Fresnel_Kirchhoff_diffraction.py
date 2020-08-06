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
    def __init__(self,Ein,Xin,Yin,Xout,Yout,L,lam=1030*10**-6):
        self.Ein=Ein
        self.Xin=Xin
        self.Yin=Yin
        self.Xout=Xout
        self.Yout=Yout
        self.L=L
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
    k=None
    
    def propagate(self):
        Xin=self.Xin #assumed to be same as Yin
        lam=self.lam
        dx=Xin[1]-Xin[0]
        dy=dx
        Nx=len(Xin)
        XinM=Xin[:,None]*np.ones(Nx)
        YinM=np.ones(Nx)[:,None]*Xin
        self.Eout=np.array([[-1j/lam*dx*dy*np.sum(self.Ein*self.R_part(XinM,YinM,x,y,self.L)*
                                                    self.Cos_part(XinM,YinM,x,y,self.L)) 
                             for y in self.Yout] for x in self.Xout])
    
    
    def R_part(self,xin,yin,xout,yout,z):
        """exp(i*k*R)/R"""
        R=((xin-xout)**2+(yin-yout)**2+z**2)**0.5
        # print(R)
        return np.exp(1j*self.k*R)/R
    
    def Cos_part(self,xin,yin,xout,yout,z):
        """(1+cos(Rz))/2"""
        # Rvec=np.array([xin-xout,yin-yout,z])
        # Rvec=Rvec/np.sum(Rvec**2)**0.5 #normalization
        # Zv=np.array([0,0,1])
        # Cos=np.dot(Rvec,Zv)
        
        Cos=z/((xin-xout)**2+(yin-yout)**2+z**2)**0.5
        
        # return (1+Cos)/2
        return Cos
        
        
        