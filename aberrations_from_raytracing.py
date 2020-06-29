"""
aberrations using ray tracing
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)

from classes.ray_tracing_class import ray_trace
from classes.ray_class import (ray_class, gaussian_angular_distribution, 
                               gaussian_spaceYangle_distribution,angular_conus,beam2aberrations)
from classes.aberration_functions import remove_low_aberrations
import numpy as np
Pi=np.pi
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy import interpolate

M_conf=[{'Mirror_type' : 'spherical', 'p' : 750, 'q' : 2, 'alfa' : 3, 'R' : 1500, 'type' : 'refocusing', 'size' : [200,100]}]

#generate rays
Rin=1.5 #radius in mm
Nrays=200 # number of rays
rays,X,Y=beam2aberrations(Rin,Nrays)
#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
trace.propagate()
trace.show_rays()
trace.show_plane()
Phase=trace.aberrations()[0].reshape((len(X),len(Y)))
Rout=min([np.abs(X).max(),np.abs(Y).max()])

P=[]

for i in range(len(X)):
    for j in range(len(Y)):
        P+=[[X[i],Y[j]]]
        if X[i]**2+Y[j]**2 > Rout**2:
            Phase[i][j]=0
            
            
FitPhase=interpolate.CloughTocher2DInterpolator(P,Phase.reshape(len(X)*len(Y)))
            
# plt.pcolormesh(X,Y,Phase)
# plt.show()

Xout=np.linspace(-Rout,Rout,int((len(X)*len(Y))**0.5))
# print(X,Xout)
Yout=np.linspace(-Rout,Rout,int((len(X)*len(Y))**0.5))
PhaseF=np.array([[FitPhase(x,y) for y in Yout] for x in Xout])


Phase1,tilt,tip,defoc=remove_low_aberrations(PhaseF,Xout,Yout,Rout)

# for i in range(len(Xout)):
#     for j in range(len(Yout)):
#         if Xout[i]**2+Yout[j]**2 > Rout**2:
#             Phase1[i][j]=0

print(tilt,tip,defoc)
plt.pcolormesh(X,Y,Phase1)
plt.show()
print(Phase.max()-Phase.min())
print(Phase1.max()-Phase1.min())