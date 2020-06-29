"""
raytracing for misalignment
@author: Slawa

all dimensions are in mm
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)

from ray_tracing_class import ray_trace
from ray_class import (ray_class, gaussian_angular_distribution, 
                               gaussian_spaceYangle_distribution,angular_conus,
                               gaussian_space_distribution)
import numpy as np
Pi=np.pi
import warnings
warnings.filterwarnings("ignore")
from tkinter import filedialog
import matplotlib.pyplot as plt


#nominal parameters (to find ideal specs); p - input arm mm; q - output arm mm; alfa - angle of icnidence relative to mirror normal
M_conf=[{'Mirror_type' : 'paraboloid', 'p' : 400, 'q' : 500, 'alfa' : 85, 'type' : 'focusing', 'size' : [200,100]}
       ]

# M_conf=[{'Mirror_type' : 'toroid', 'p' : 1000, 'q' : 1000, 'alfa' : 85, 'type' : 'focusing', 'size' : [300,150]}
#         ]

#rays for simulations. For point spred function simulations. all rays from same point with gaussian angular distribution
W=2.5*0.9/2#in mm; radius on the e**-2 intensity level 
Nrays=400 # number of rays (could be slightly reduced to acceletare computations)
#ganerate rays for simulations
rays=gaussian_space_distribution(W,Nrays)


#%%
#========== angular tangential misalignment (one point example)==============
#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)
Mtilt=150*10**-6 #radians mirror tilt in the reflection plane
trace.mirrors[0].tiltx(Mtilt)
# trace.mirrors[0].q+=-35
trace.config(skip_align=True) #find new output plane

trace.propagate()
# trace.show_rays()
trace.show_plane()
print(trace.Mray.X)

#%%
#========== angular saggital misalignment (one point example)==============

p=M_conf[0]['p']
Mtilt=150*10**-6 #radians mirror tilt in the reflection plane
X0=[p-p*np.cos(Mtilt),0,-p*np.sin(Mtilt)]
A0=[np.cos(Mtilt),0,np.sin(Mtilt)]
rays=gaussian_space_distribution(W,Nrays,X0,A0)
# for r in rays:
#     r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf,MrayX=X0,MrayA=A0)

# trace.mirrors[0].tiltx(Mtilt)
# trace.mirrors[0].q+=-35
# trace.config(skip_align=True) #find new output plane

trace.propagate()
trace.show_rays()
trace.show_plane()
print(trace.Mray.X)

#%%
#============shift perpendicular to the input beam (one point example)==========
#in the input plane (the plane between the k of the input beam and the normal to the surface)

#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)

Lshift=0.08 #mm amount of shift
m=trace.mirrors[0]
m.shift([0,Lshift,0])
#find new output plane

q2=trace.find_focus(collimated=True)
#apply new output arm length
# m.q=q2
print(q2)
# m.q+=0.5

trace.config(skip_align=True) #find new output plane
trace.propagate()
# trace.show_rays()
trace.show_plane()

#%%
#=====================shift vertically (one point example)==================


#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)

Lshift=50 #mm amount of shift
m=trace.mirrors[0]
m.shift([0,0,Lshift])
#find new output plane

q2=trace.find_focus(collimated=True)
#apply new output arm length
m.q=q2
print(q2)
# m.q+=0.5

trace.config(skip_align=True) #find new output plane
trace.propagate()
# trace.show_rays()
trace.show_plane()


#%%






#=========================scan example==============================

#=================angular tangential misalignment (scan example with many points)===============
#scans for other misalignments can be done in the same way

#angles to scan
Mtilt=np.linspace(-5*10**-3,5*10**-3,11)
dmt=Mtilt[1]-Mtilt[0]

Output=[] #list for the output data

#prepare rays and mirrors
for r in rays:
    r.reset()
trace=ray_trace(rays,mirrors_config=M_conf)
trace.mirrors[0].tiltx(Mtilt[0]-dmt)

for i in range(len(Mtilt)):
    for r in rays:
        r.reset()
    trace.mirrors[0].tiltx(dmt)
    trace.config(skip_align=True)
    trace.propagate()
    Output.append(trace.rms())

Rrms=[ot[2]*10**3 for ot in Output]
plt.plot(Mtilt*10**3,Rrms)
plt.xlabel('tilt (mrad)')
plt.ylabel('PSF (um)')
plt.show()


#%%
#=================angular saggital misalignment (scan example with many points)===============
#scans for other misalignments can be done in the same way

#angles to scan
Mtilt=np.linspace(-5*10**-3,5*10**-3,11)
dmt=Mtilt[1]-Mtilt[0]

Output=[] #list for the output data

#prepare rays and mirrors
p=M_conf[0]['p']



for i in range(len(Mtilt)):
    
    Mt=Mtilt[i]
    X0=[p-p*np.cos(Mt),0,-p*np.sin(Mt)]
    A0=[np.cos(Mt),0,np.sin(Mt)]
    rays=gaussian_space_distribution(W,Nrays,X0,A0)
    trace=ray_trace(rays,mirrors_config=M_conf,MrayX=X0,MrayA=A0)
    # trace.config(skip_align=True)
    trace.propagate()
    Output.append(trace.rms())

Rrms=[ot[2]*10**3 for ot in Output]
plt.plot(Mtilt*10**3,Rrms)
plt.xlabel('tilt (mrad)')
plt.ylabel('PSF (um)')
plt.show()