"""
raytracing for misalignment
@author: Slawa

all dimensions are in mm
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)
sys.path.append(Path+'//classes')

from ray_tracing_class import ray_trace
from ray_class import (ray_class, gaussian_angular_distribution, 
                               gaussian_spaceYangle_distribution,angular_conus)
import numpy as np
Pi=np.pi
import warnings
warnings.filterwarnings("ignore")
from tkinter import filedialog
import matplotlib.pyplot as plt

#nominal parameters (to find ideal specs); p - input arm mm; q - output arm mm; alfa - angle of icnidence relative to mirror normal
M_conf=[{'Mirror_type' : 'plane_ellips_z', 'p' : 2450, 'q' : 700, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
        ]

# M_conf=[{'Mirror_type' : 'ellipsoid', 'p' : 1000, 'q' : 1000, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
#        ]

#rays for simulations. For point spred function simulations. all rays from same point with gaussian angular distribution
AngularWidth=2*10**-3 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
Nrays=400 # number of rays (could be slightly reduced to acceletare computations)
#ganerate rays for simulations
rays=gaussian_angular_distribution(AngularWidth,Nrays)

#%%
#========== angular misalignment (one point example)==============
#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)
Mtilt=50*10**-4 #radians mirror tilt in the reflection plane
trace.mirrors[0].tiltx(Mtilt)
p=M_conf[0]['p']
q=M_conf[0]['q']
alfa=M_conf[0]['alfa']/180*Pi
f=1/(1/p+1/q)
f1=f/np.cos(alfa+Mtilt)*np.cos(alfa)
print('focus shift ',f1-f)
# q1=1/(1/f1-1/p)

q1=trace.find_focus()


trace.mirrors[0].q=q1
# trace.mirrors[0].q+=3
trace.config(skip_align=True) #find new output plane

trace.propagate()
# trace.show_rays()
trace.show_plane()


#%%







#=========================scan example==============================

#=================angular tangential misalignment (scan example with many points)===============

#angles to scan
Mtilt0=np.linspace(-15*10**-3,15*10**-3,11)

Output=[] #list for the output data
p=M_conf[0]['p']
q=M_conf[0]['q']
f=1/(1/p+1/q)
alfa=M_conf[0]['alfa']/180*Pi
#prepare rays and mirrors
for r in rays:
    r.reset()

for i in range(len(Mtilt0)):
    for r in rays:
        r.reset()
    Mtilt=Mtilt0[i]
    trace=ray_trace(rays,mirrors_config=M_conf)
    trace.mirrors[0].tiltx(Mtilt)
    # f1=f/np.cos(alfa+Mtilt)*np.cos(alfa)
    # q1=1/(1/f1-1/p)
    q1=trace.find_focus()
    trace.mirrors[0].q=q1
    trace.config(skip_align=True)
    trace.propagate()
    Output.append(trace.rms())

Rrms=np.array([ot[0]*10**3 for ot in Output])
plt.plot(Mtilt0*10**3,Rrms*4)
plt.xlabel('tilt (mrad)')
plt.ylabel('PSF: beam diameter (um)')
plt.show()

Rrms1=np.concatenate((Rrms[-1:-6:-1],Rrms[5:]))
Rrms=Rrms1

file=filedialog.asksaveasfilename()
if file != '':
    np.savetxt(file+'.dat', np.concatenate((np.array(Mtilt0*10**3).reshape((-1,1)), np.array(Rrms*4).reshape((-1,1))),axis=1),
        delimiter='\t',comments='')