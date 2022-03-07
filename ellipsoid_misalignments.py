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
from tkinter import filedialog


#nominal parameters (to find ideal specs); p - input arm mm; q - output arm mm; alfa - angle of icnidence relative to mirror normal
M_conf=[{'Mirror_type' : 'ellipsoid', 'p' : 2250, 'q' : 900, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
        ]

# M_conf=[{'Mirror_type' : 'ellipsoid', 'p' : 1000, 'q' : 1000, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
#        ]

#rays for simulations. For point spred function simulations. all rays from same point with gaussian angular distribution
AngularWidth=2.5*10**-3*1.8/2 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
Nrays=400 # number of rays (could be slightly reduced to acceletare computations)
#ganerate rays for simulations
rays=gaussian_angular_distribution(AngularWidth,Nrays)


#%%
#========== angular tangential misalignment (one point example)==============
#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)
Mtilt=25*10**-6 #radians mirror tilt in the reflection plane
trace.mirrors[0].tiltx(Mtilt)
# trace.mirrors[0].q+=5
trace.config(skip_align=True) #find new output plane

trace.propagate()
# trace.show_rays()
trace.show_plane()
# trace.export_plane()

#%%
#========== angular saggital misalignment (one point example)==============

p=M_conf[0]['p']
Mtilt=25*10**-6 #radians mirror tilt in the reflection plane
X0=[p-p*np.cos(Mtilt),0,-p*np.sin(Mtilt)]
A0=[np.cos(Mtilt),0,np.sin(Mtilt)]
rays=gaussian_angular_distribution(AngularWidth,Nrays,X0,Mtilt)
# for r in rays:
#     r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf,MrayX=X0,MrayA=A0)

# trace.mirrors[0].tiltx(Mtilt)
# trace.mirrors[0].q+=-5
# trace.config(skip_align=True) #find new output plane

trace.propagate()
# trace.show_rays()
trace.show_plane()
# print(trace.Mray.X)

#%%
#=====================shift vertically (one point example)==================
#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)

Lshift=0.055 #mm
m=trace.mirrors[0]
m.shift([0,0,Lshift])
print('\n should be about equivalent to the vertical tillt of' , Lshift/m.r, 'rad \n')

#find new output plane (and the length of the output arm)

p=M_conf[0]['p']
q=M_conf[0]['q']
f=1/(1/p+1/q) #assuming refocusing
#find new input arm length
XYZ=m.surf_coordinates([[0],[0]],given = 'yz') #intersection coodinates with the surface
p2=0
X0=m.XM0[0]
#find the correct intersection
for a in XYZ:
    if np.abs(a[0]-X0) < Lshift*4+10**-6:
        p2=float(a[0][0])
# q2=1/(1/f-1/p2)
q2=q+(p-p2)
#apply new output arm length
m.q=q2
m.q+=0

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

Lshift=0.05 #mm amount of shift
m=trace.mirrors[0]
m.shift([0,Lshift,0])
#find new output plane

p=M_conf[0]['p']
q=M_conf[0]['q']
f=1/(1/p+1/q) #assuming refocusing
# #find new input arm length
XYZ=m.surf_coordinates([[0],[0]],given = 'yz') #intersection coodinates with the surface
p2=0
X0=m.XM0[0]
# #find the correct intersection
for a in XYZ:
    if a[0]-X0 < Lshift*4+10**-6:
        p2=float(a[0][0])
# #new angle of incidence
# # alfa2=np.arctan(m.X0[1]/(p2-m.X0[0]))
# # f2=f*np.cos(alfa2)/np.cos(M_conf[0]['alfa']/180*Pi) #new focal length
# #new output arm length
# q2=1/(1/f-1/p2)

q2=trace.find_focus()
# q2=q+(p-p2)
#apply new output arm length
m.q=q2
m.q+=0

trace.config(skip_align=True) #find new output plane
trace.propagate()
# trace.show_rays()
trace.show_plane()

#%%
# ====================shift alog the input beam (one point example)================
#reset rays to original conditions
for r in rays:
    r.reset()
#create the trace class and get ideal parameters
trace=ray_trace(rays,mirrors_config=M_conf)

Lshift=30 #mm amount of shift
trace.mirrors[0].shift([Lshift,0,0])
#find new output arm length
p=M_conf[0]['p']
q=M_conf[0]['q']
f=1/(1/p+1/q) #assuming refocusing
p2=p+Lshift
q2=1/(1/f-1/p2)
print('new out arm: ', q2,' total length change: ', p+q-(p2+q2))
# q2=q-Lshift
#apply new output arm length
trace.mirrors[0].q=q2

trace.config(skip_align=True) #find new output plane
trace.propagate()
# trace.show_rays()
trace.show_plane()

#%%







#=========================scan example==============================

#=================angular tangential  misalignment (scan example with many points)===============

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
plt.xlabel('tilt ($\mu$rad)')
plt.ylabel('PSF: beam diameter ($\mu$m)')
plt.show()

# file=filedialog.asksaveasfilename()
# if file != '':
#     np.savetxt(file+'.dat', np.concatenate((np.array(Mtilt*10**3).reshape((-1,1)), np.array(Rrms*4).reshape((-1,1))),axis=1),
#         delimiter='\t',comments='')

#%%
#=================angular saggital misalignment (scan example with many points)===============

#angles to scan
Mtilt=np.linspace(-100*10**-6,100*10**-6,11)

Output=[] #list for the output data

#prepare rays and mirrors
p=M_conf[0]['p']



for i in range(len(Mtilt)):
    
    Mt=Mtilt[i]
    X0=[p-p*np.cos(Mt),0,-p*np.sin(Mt)]
    A0=[np.cos(Mt),0,np.sin(Mt)]
    rays=gaussian_angular_distribution(AngularWidth,Nrays,X0,Mt)
    trace=ray_trace(rays,mirrors_config=M_conf,MrayX=X0,MrayA=A0)
    # trace.config(skip_align=True)
    trace.propagate()
    Output.append(trace.rms())

Rrms=np.array([ot[2]*10**3 for ot in Output])
plt.plot(Mtilt*10**6,Rrms*4)
plt.xlabel('tip ($\mu$rad)') #rather pitch
plt.ylabel('PSF: beam diameter ($\mu$m)')
plt.show()

# file=filedialog.asksaveasfilename()
# if file != '':
#     np.savetxt(file+'.dat', np.concatenate((np.array(Mtilt*10**3).reshape((-1,1)), np.array(Rrms*4).reshape((-1,1))),axis=1),
#         delimiter='\t',comments='')

#%%
#=================vertical shift misalignment (scan example with many points)===============

#shifts to scan
Ls=np.linspace(-0.2,0.2,11)
# dl=Ls[1]-Ls[0]

Output=[] #list for the output data
p=M_conf[0]['p']
q=M_conf[0]['q']
f=1/(1/p+1/q) #assuming refocusing

for i in range(len(Ls)):
    for r in rays:
        r.reset()
    trace=ray_trace(rays,mirrors_config=M_conf)
    Lshift=Ls[i]
    m=trace.mirrors[0]
    m.shift([0,0,Lshift])

    #find new input arm length
    XYZ=m.surf_coordinates([[0],[0]],given = 'yz') #intersection coodinates with the surface
    p2=0
    X0=m.XM0[0]
    #find the correct intersection
    for a in XYZ:
        if np.abs(a[0]-X0) > 0:
            p2=float(a[0][0])
    q2=1/(1/f-1/p2)
    #apply new output arm length
    m.q=q2
    
    trace.config(skip_align=True)
    trace.propagate()
    Output.append(trace.rms())

Rrms=np.array([ot[2]*10**3 for ot in Output])
plt.plot(Ls*10**3,Rrms*4)
plt.xlabel('vertival shift ($\mu$m)')
plt.ylabel('PSF RMS diameter ($\mu$m)')
plt.show()



