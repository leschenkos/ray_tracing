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
                               gaussian_spaceYangle_distribution,angular_conus,
                               gaussian_space_distribution)
import numpy as np
Pi=np.pi
import warnings
warnings.filterwarnings("ignore")
from tkinter import filedialog
import matplotlib.pyplot as plt


#nominal parameters (to find ideal specs); p - input arm mm; q - output arm mm; alfa - angle of icnidence relative to mirror normal
M_conf=[{'Mirror_type' : 'paraboloid', 'p' : 900, 'q' : 150, 'alfa' : 85, 'type' : 'collimating', 'size' : [200,100]}
        ,{'Mirror_type' : 'paraboloid', 'p' : 150, 'q' : 500, 'alfa' : 85, 'type' : 'focusing', 'size' : [200,100]}
       ]

# M_conf=[{'Mirror_type' : 'toroid', 'p' : 1000, 'q' : 1000, 'alfa' : 85, 'type' : 'focusing', 'size' : [300,150]}
#         ]

# #rays for simulations. For point spred function simulations. all rays from same point with gaussian angular distribution
# W=0.5#in mm; radius on the e**-2 intensity level 
# Nrays=400 # number of rays (could be slightly reduced to acceletare computations)
# #ganerate rays for simulations
# rays=gaussian_space_distribution(W,Nrays)

#rays for simulations. For point spred function simulations. all rays from same point with gaussian angular distribution
AngularWidth=2.5*10**-3/2 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
Nrays=400 # number of rays (could be slightly reduced to acceletare computations)
#ganerate rays for simulations
rays=gaussian_angular_distribution(AngularWidth,Nrays)


#%%
#========== first mirror -> angular tangential misalignment (one point example)==============
#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)
Mtilt=250*10**-6*1 #radians mirror tilt in the reflection plane
trace.mirrors[0].tiltx(Mtilt)
# trace.mirrors[0].q+=-35
trace.config(skip_align=True) #find new output plane

trace.propagate()
# trace.show_rays()
trace.show_plane()
print(trace.Mray.X)

#%%
#========== second mirror -> angular tangential misalignment (one point example)==============
#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)
Mtilt=150*10**-6*1 #radians mirror tilt in the reflection plane
trace.mirrors[1].tiltx(Mtilt)
# trace.mirrors[0].q+=-35
trace.config(skip_align=True,skipmirror=1) #find new output plane

trace.propagate()
# trace.show_rays()
trace.show_plane()
print(trace.Mray.X)

#%%
#========== first mirror -> shift alog the input beam (one point example)====
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)

Lshift=3 #mm amount of shift
trace.mirrors[0].shift([Lshift,0,0])
# print(trace.mirrors[1].XM0)
# trace.config(skip_align=True)
# # print(trace.mirrors[1].XM0)
# q2=trace.find_focus()
# print(q2)
# # trace.mirrors[1].q=q2
# trace.mirrors[1].q+=0.32

p=M_conf[0]['p']
q=M_conf[1]['q']
f=1/(1/p+1/q) #assuming refocusing
p2=p+Lshift
q2=1/(1/f-1/p2)
# trace.mirrors[1].q=2*q-q2
# trace.mirrors[1].q+=7
print(trace.mirrors[1].q)

trace.config(skip_align=True) #find new output plane
trace.propagate()
# trace.show_rays()
trace.show_plane()


#%%
#============shift perpendicular to the input beam (one point example)==========
#in the input plane (the plane between the k of the input beam and the normal to the surface)

#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)

Lshift=0.2 #mm amount of shift
m=trace.mirrors[0]
m.shift([0,Lshift,0])
#find new output plane
trace.config(skip_align=True) 
q2=trace.find_focus()
print(q2)
#apply new output arm length
trace.mirrors[1].q=q2
trace.mirrors[1].q+=0

trace.config(skip_align=True) #find new output plane
trace.propagate()
# trace.show_rays()
trace.show_plane()

#%%











#=========================scan example==============================

#=================angular tangential misalignment (scan example with many points)===============
#scans for other misalignments can be done in the same way

#angles to scan
Mtilt=np.linspace(-100*10**-6,100*10**-6,11)
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

Rrms=np.array([ot[2]*10**3 for ot in Output])
plt.plot(Mtilt*10**6,Rrms*4)
plt.xlabel('tilt Par_1 ($/mu$rad)')
plt.ylabel('PSF ($/mu$m)')
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