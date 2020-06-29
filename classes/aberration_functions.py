"""
functions to work with aberrations
https://en.wikipedia.org/wiki/Zernike_polynomials
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)

import numpy as np
Pi=np.pi
import matplotlib.pyplot as plt

def Zer_x(x,y):
    """Zernike tilt polinom"""
    r=(x**2+y**2)**0.5
    if r <= 1:
        return x
    else:
        return 0

def Zer_y(x,y):
    """Zernike tip polinom"""
    r=(x**2+y**2)**0.5
    if r <= 1:
        return y
    else:
        return 0

def Zer_def(x,y):
    """Zernike defocus polinom"""
    r=(x**2+y**2)**0.5
    if r<=1:
        return 2*r**2-1
    else:
        return 0
    
def remove_low_aberrations(Phase0,X,Y,R):
    """remove tilt, tip and defocus"""
    Zero=mean_round(Phase0,X,Y,R)
    Phase=zero_out(Phase0-Zero,X,Y,R)
    # R=max(np.abs(X).max(),np.abs(Y).max()) #radius
    Zx=np.array([[Zer_x(x/R,y/R) for y in Y] for x in X])
    Zy=np.array([[Zer_y(x/R,y/R) for y in Y] for x in X])
    Zf=np.array([[Zer_def(x/R,y/R) for y in Y] for x in X])
    Zx=zero_out(Zx-mean_round(Zx,X,Y,R),X,Y,R)
    Zy=zero_out(Zy-mean_round(Zy,X,Y,R),X,Y,R)
    Zf=zero_out(Zf-mean_round(Zf,X,Y,R),X,Y,R)
    tilt=np.sum(Phase*Zx)/np.sum(Zx*Zx)
    tip=np.sum(Phase*Zy)/np.sum(Zy*Zy)
    defoc=np.sum(Phase*Zf)/np.sum(Zf*Zf)
    Phase+=-tilt*Zx-tip*Zy-defoc*Zf
    # print(mean_round(Phase,X,Y,R))
    # print(np.sum(Phase*Zf)/np.sum(Zf*Zf))
    return Phase,tilt,tip,defoc

def mean_round(Phase,X,Y,R):
    """mean value in circle
    phase need to be NxN np.array"""
    Pr=[]
    for i in range(len(X)):
        for j in range(len(Y)):
            if (X[i]**2+Y[j]**2) <= R**2:
                Pr+=[Phase[i][j]]
    return np.array(Pr).mean()

    # N=len(Phase)
    # R=(N-1)/2
    # Pr=[]
    # X=np.linspace(0,N-1,N)-R
    # Y=X
    # for i in range(N):
    #     for j in range(N):
    #         if (X[i]**2+Y[j]**2) <= R**2:
    #            Pr+=[Phase[i][j]] 
    # return np.array(Pr).mean()
    
def zero_out(Phase,X,Y,R):
    """zero phase outside R"""
    for i in range(len(X)):
        for j in range(len(Y)):
            if X[i]**2+Y[j]**2 > R**2 + 10**-13:
                Phase[i][j]=0
    return Phase
    
    
    