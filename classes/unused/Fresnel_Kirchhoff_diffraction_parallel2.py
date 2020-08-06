"""
Fresnel-Kirchhoff diffraction
all distances in mm

to do:
parallel computing
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
        Xout=self.Xout
        
        #to do:
            #split into Ncore subtasks
        # Nc=p._processes
        Nc=11
        Nout=len(self.Xout)
        nn=int(Nout/Nc) #number of points per process
        # Nn=np.ones(Nc)
        # Nn[:-2]=nn
        # Nn[-1]=Nout-nn*(Nc-1)
        Nn=[list(range(nn*i,nn*(i+1))) for i in range(Nc-1)]
        Nn+= [list(range(nn*(Nc-1),Nout))]
        # print(Nn)
        E=p.starmap(prop_x,[[Ein,XinM,YinM,Xout,Yout,L,lam,dx,nn] for nn in Nn])
        E2=[]
        for e in E:
            E2.extend(e)
        self.Eout=np.array(E2)
        
    
def prop_x(Ein,Xin,Yin,Xout,Yout,L,lam,dx,Nx):
    """for parallelzation
    propagates 1 raw of data at fixed x position"""
    dy=dx
    Xout1=Xout[Nx]
    Out=[[-1j/lam*dx*dy*np.sum(Ein*R_part(Xin,Yin,x,y,L,2*Pi/lam)*
                                              Cos_part(Xin,Yin,x,y,L))
                  for y in Yout]
         for x in Xout1]
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
    
    
# F=FK_diffraction(1,np.ones(100),np.ones(100),np.ones(100),np.ones(100),1,1)
# F.propagate()