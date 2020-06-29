"""
tilted mirror aberrations

@author: Slawa
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)

from ray_tracing_class import ray_trace
from ray_class import (ray_class, gaussian_angular_distribution, 
                               gaussian_spaceYangle_distribution,angular_conus,beam2aberrations)
from aberration_functions import remove_low_aberrations
import numpy as np
Pi=np.pi
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import interpolate

def mirror_aberration(R,alfa,X,Y,lam=1030*10**-6):
    """returns phase distortions of a spherical mirror with radius R
    alfa tilt angle in degree
    X, Y is a space mask to compute aberrations (they are expected to be symmetrical around 0, and the mirror is assumed to be round)
    lam wavelength in mm"""
    
    M_conf=[{'Mirror_type' : 'spherical', 'p' : 500, 'q' : 2, 'alfa' : alfa, 'R' : R, 'type' : 'refocusing', 'size' : [200,100]}]
    
    rays=beam2aberrations(1,1,X,Y)
    trace=ray_trace(rays,mirrors_config=M_conf)
    trace.propagate()
    Phase=trace.aberrations(lam)[0].reshape((len(X),len(Y)))
    plt.pcolormesh(X,Y,Phase)
    plt.show()
    # q=M_conf[0]['q']
    # plt.pcolormesh(X,Y,Phase-defocus(R-q,X,Y,lam))
    # plt.show()
    Rout=min([np.abs(X).max(),np.abs(Y).max()]) #radius of the mirror (aberrations are fitted in the circle)
    # print(Rout)
    for i in range(len(X)):
        for j in range(len(Y)):
            if X[i]**2+Y[j]**2 > Rout**2 + 10**-13: #10**-13 to get rid of troubles right on the border of the radius
                Phase[i][j]=0
    # plt.pcolormesh(X,Y,Phase)
    # plt.show()
    Phase1,tilt,tip,defoc=remove_low_aberrations(Phase,X,Y,Rout)
    return Phase1

def defocus(R,X,Y,lam=1030*10**-6):
    k=2*Pi/lam
    return np.array([[-k/R*(x**2+y**2) for y in Y] for x in X])

def telescope_aberration(R1,alfa1,R2,alfa2,X,Y,lam=1030*10**-6):
    """returns phase distortions of a spherical mirror with radius R; >0 focusing, <0 defocusing
    alfa tilt angle in degree; 1 : first mirror; 2 : second mirror
    X, Y is a space mask to compute aberrations (they are expected to be symmetrical around 0, and the mirror is assumed to be round)
    lam wavelength in mm"""
    
    f1=R1/2
    f2=R2/2
    if f1+f2 < 0:
        print('negative f1+f2; wrong telescope parameters')
    M_conf=[{'Mirror_type' : 'spherical', 'p' : np.abs(f1), 'q' : (f1+f2)/2, 'alfa' : alfa1, 'R' : R1, 'type' : 'refocusing', 'size' : [200,100]},
            {'Mirror_type' : 'spherical', 'p' : (f1+f2)/2, 'q' : np.abs(f2), 'alfa' : alfa2, 'R' : R2, 'type' : 'refocusing', 'size' : [200,100]}]
    
    rays=beam2aberrations(1,1,X,Y)
    trace=ray_trace(rays,mirrors_config=M_conf)
    trace.propagate()
    Phase=trace.aberrations(lam)[0].reshape((len(X),len(Y)))
    # plt.pcolormesh(X,Y,Phase)
    # plt.show()
    # q=M_conf[0]['q']
    # plt.pcolormesh(X,Y,Phase-defocus(R-q,X,Y,lam))
    # plt.show()
    Rout=min([np.abs(X).max(),np.abs(Y).max()]) #radius of the mirror (aberrations are fitted in the circle)
    # print(Rout)
    for i in range(len(X)):
        for j in range(len(Y)):
            if X[i]**2+Y[j]**2 > Rout**2 + 10**-13: #10**-13 to get rid of troubles right on the border of the radius
                Phase[i][j]=0
    # plt.pcolormesh(X,Y,Phase)
    # plt.show()
    Phase1,tilt,tip,defoc=remove_low_aberrations(Phase,X,Y,Rout)
    return Phase1

# X=np.linspace(-4,4,81)
# Y=X
# Ph=mirror_aberration(1500,2,X,Y)

# plt.pcolormesh(X,Y,Ph)
# plt.show()

# X=np.linspace(-4,4,81)
# Y=X
# Ph=telescope_aberration(1500,2,-750,-2,X,Y)

# plt.pcolormesh(X,Y,Ph)
# plt.show()
    
    
    
    