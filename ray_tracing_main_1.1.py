"""
ray tracing v 1.3 // 08.05.20
@author: Slawa

all dimensions are in mm

it a 3D code based on geometrical optics approximations
(which should be perfectly valid for most of the practical cases of XUV propagation,
 since typical focal diameter is much larger than a wavelength)
each ray propagates straight till a surface
on each surface angle of reflection equal to the angle of incidence (relative to the normal of the surface at the intersection point)

the mirrors config file is a list of dictionaries (see examples in the file)

{'Mirror_type' : 'toroid', 'p' : 1000, 'q' : 1000, 'alfa' : 80, 'type' : 'refocusing', 'size' : [200,100]}
right know: Mirror_type can be: toroid, ellipsoid, spherical and plane
p - length of the input arm
q - length of the output arm
alfa - angle of incidence in degrees (~80 is a grazing incidence)
positive angle means reflection to the left relative to the incoming light; negative to the right
the mirror parameters (i.e. radius) are computed automatically according to the above specs and the type
type: for toroid: 
    refocusing (refocuses a source from p to q), 
    collimating (collimates a point source at p), 
    focusing (focuses a collimated beam at q), and 
    fixed (radii R and r, R>=r, need to be specified)
        fixed mode can be used to study the alignment accuracy requirements
type: for ellipsoid:
    refocusing (refocuses a source from p to q),
    fixed (radii R and r, R>=r, need to be specified)
spherical and plane have only one type: 'fixed' (which doesn't have to be specified)
size is the size of the mirror [width,hight]

presently there are three options to generate the input set of rays
1. angular_conus : conical set of rays starting from the same point
can be used for rough and very fast estimations
2. gaussian_angular_distribution : a gaussian like angular distribution of the rays starting from the same point
requires more rays and more computation time (still 10s of seconds). 
can be used for point spread function evaluation
3. gaussian_spaceYangle_distribution : generates gaussian like spatial and angular distribution
each starting point has a gaussian like angular distribution (identical to gaussian_angular_distribution)

any other ray distribution can be easily implemented (contact me if you have troubles)

the final data representation is very simple right know (contact me leshchenko.1@osu.edu 
                                    if you have troubles with getting out any extra data)


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


#%%
"""quick check of a multi mirror configuration
provides a rough idea of what to expect (gives the outer borders of the beam)
"""

M_conf=[{'Mirror_type' : 'ellipsoid', 'p' : 1500, 'q' : 500, 'alfa' : 80, 'type' : 'refocusing', 'size' : [200,100]}
        ,{'Mirror_type' : 'toroid', 'p' : 20, 'q' : 400, 'alfa' : -80, 'type' : 'collimating', 'size' : [200,100]}
        ,{'Mirror_type' : 'toroid', 'p' : 300, 'q' : 600, 'alfa' : 80, 'type' : 'focusing', 'size' : [200,100]}
        ]

#generate rays
ConAng=1*10**-3 #in rad; angle of the conus
Nrays=100 # number of rays
rays=angular_conus(ConAng,Nrays)
#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
trace.propagate()
trace.show_rays()
trace.show_plane()

#%%
"""more fancy ray set (it is an example similar to Opt. Express 21, 13040-13051 (2013), the overall size is comparable, though the focal distribution looks different)"""

M_conf=[{'Mirror_type' : 'toroid', 'p' : 1500, 'q' : 150, 'alfa' : 80, 'type' : 'collimating', 'size' : [200,100]}
        ,{'Mirror_type' : 'toroid', 'p' : 150, 'q' : 150, 'alfa' : -80, 'type' : 'focusing', 'size' : [200,100]}
        ,{'Mirror_type' : 'toroid', 'p' : 660, 'q' : 600, 'alfa' : 80, 'type' : 'refocusing', 'size' : [200,100]}
        ]

#generate rays
AngularWidth=1*10**-3 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
Nrays=1000 # number of rays
rays=gaussian_angular_distribution(AngularWidth,Nrays)
#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
trace.propagate()
trace.show_rays()
trace.show_plane()

#%%
"""more fancy ray set
point spread function assuming gaussian angular distribution
the toroid parameters are similar to the FRED simulations in the Greg's thesis"""
M_conf=[{'Mirror_type' : 'toroid', 'p' : 1500, 'q' : 500, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
        ]

#generate rays
AngularWidth=2*10**-3 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
Nrays=1000 # number of rays
rays=gaussian_angular_distribution(AngularWidth,Nrays)
#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
trace.propagate()
trace.show_rays()
trace.show_plane()



#%%
"""most fancy rays set
example of refocusing a kind of Gaussian beam (with a set of correspondingly distributed rays)
the number of rays is quite small (to make it fast), so the picture looks like a snowflake,
 but it is already enough to get a filling of an expected size"""

AngularWidth=2*10**-3 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
Wradius=0.1 #starting beam radius on the e**-2 intensity level 
Nrays=4000 # number of rays
rays=gaussian_spaceYangle_distribution(Wradius,int(np.sqrt(Nrays)),AngularWidth,int(np.sqrt(Nrays)),Nsub0=6)

M_conf=[{'Mirror_type' : 'ellipsoid', 'p' : 1000, 'q' : 300, 'alfa' : 80, 'type' : 'refocusing', 'size' : [200,100]}
        ,{'Mirror_type' : 'toroid', 'p' : 1000, 'q' : 1000, 'alfa' : -80, 'type' : 'refocusing', 'size' : [200,100]}
        ]

trace=ray_trace(rays,mirrors_config=M_conf)
trace.propagate()
trace.show_rays()
trace.show_plane()

#%%
"""point spred function of a monochromator with paraboliods
(shows numerical accuracy since in theory it should end up with a perfect point)
"""

M_conf=[{'Mirror_type' : 'paraboloid', 'p' : 1000, 'q' : 150, 'alfa' : 85, 'type' : 'collimating', 'size' : [200,100]}
        ,{'Mirror_type' : 'paraboloid', 'p' : 150, 'q' : 500, 'alfa' : -85, 'type' : 'focusing', 'size' : [200,100]}
        ]

#generate rays
AngularWidth=1*10**-3 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
Nrays=1000 # number of rays
rays=gaussian_angular_distribution(AngularWidth,Nrays)

#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
trace.propagate()
trace.show_rays()
trace.show_plane()