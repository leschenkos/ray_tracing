"""
raytracing for misalignment
@author: Slawa

all dimensions are in mm
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)
# SP=Path.split("\\")
# i=0
# while i<len(SP) and SP[i].find('python')<0:
#     i+=1
# Pypath='\\'.join(SP[:i+1])
# sys.path.append(Pypath)

from classes.ray_tracing_class import ray_trace
from classes.ray_class import (ray_class, gaussian_angular_distribution, 
                               gaussian_spaceYangle_distribution,angular_conus)
import numpy as np
Pi=np.pi
import warnings
warnings.filterwarnings("ignore")
from tkinter import filedialog
import matplotlib.pyplot as plt

"""misalignment of a toroidal mirror"""
#nominal parameters (to find ideal specs); p - input arm mm; q - output arm mm; alfa - angle of icnidence relative to mirror normal
M_conf=[{'Mirror_type' : 'toroid', 'p' : 800, 'q' : 800, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
       ]

#rays for simulations. For point spred function simulations. all rays from same point with gaussian angular distribution
AngularWidth=5*10**-3 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
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
Mtilt=-5*10**-4 #radians mirror tilt in the reflection plane
trace.mirrors[0].tiltx(Mtilt)
trace.config(skip_align=True) #find new output plane

trace.propagate()
trace.show_rays()
trace.show_plane()

#%%
# ====================shift alog the input beam (one point example)================
#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)

Lshift=-20 #mm amount of shift
trace.mirrors[0].shift([Lshift,0,0])
#find new output arm length
p=M_conf[0]['p']
f=p/2 #assuming 4f refocusing configuration
p2=p+Lshift
q2=1/(1/f-1/p2)
#apply new output arm length
trace.mirrors[0].q=q2

trace.config(skip_align=True) #find new output plane
trace.propagate()
trace.show_rays()
trace.show_plane()

#%%
#============shift perpendicular to the input beam (one point example)==========
#in the input plane (the plane between the k of the input beam and the normal to the surface)

#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)

Lshift=1 #mm amount of shift
m=trace.mirrors[0]
m.shift([0,Lshift,0])
#find new output plane

p=M_conf[0]['p']
f=p/2 #assuming 4f refocusing configuration
#find new input arm length
XYZ=m.surf_coordinates([[0],[0]],given = 'yz') #intersection coodinates with the surface
p2=0
X0=m.XM0[0]
#find the correct intersection
for a in XYZ:
    if a[0]-X0 < Lshift*4:
        p2=float(a[0][0])
#new angle of incidence
# alfa2=np.arctan(m.X0[1]/(p2-m.X0[0]))
# f2=f*np.cos(alfa2)/np.cos(M_conf[0]['alfa']/180*Pi) #new focal length
#new output arm length
q2=1/(1/f-1/p2)
#apply new output arm length
m.q=q2

trace.config(skip_align=True) #find new output plane
trace.propagate()
trace.show_rays()
trace.show_plane()

#%%
#=====================shift vertically (one point example)==================
#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)

Lshift=0.1 #mm
m=trace.mirrors[0]
m.shift([0,0,Lshift])
print('\n should be about equivalent to the vertical tillt of' , Lshift/m.r, 'rad \n')

#find new output plane (and the length of the output arm)

p=M_conf[0]['p']
f=p/2 #assuming 4f refocusing configuration
#find new input arm length
XYZ=m.surf_coordinates([[0],[0]],given = 'yz') #intersection coodinates with the surface
p2=0
X0=m.XM0[0]
#find the correct intersection
for a in XYZ:
    if np.abs(a[0]-X0) < Lshift*4+10**-6:
        p2=float(a[0][0])
q2=1/(1/f-1/p2)
#apply new output arm length
m.q=q2

trace.config(skip_align=True) #find new output plane
trace.propagate()
trace.show_rays()
trace.show_plane()

#%%
#=========================scan example==============================

#=================angular misalignment (scan example with many points)===============
#scans for other misalignments can be done in the same way

#angles to scan
Mtilt=np.linspace(-1*10**-3,1*10**-3,11)
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
plt.plot(Mtilt*10**3,Rrms*4) #RMS needs to be multiplied by 4 to get the equivalent of e**-2 diameter
plt.xlabel('tilt (mrad)')
plt.ylabel('PSF (um)')
plt.show()








