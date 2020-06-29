import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)

from classes.ray_tracing_class import ray_trace
from classes.ray_class import (ray_class, gaussian_angular_distribution, 
                               gaussian_spaceYangle_distribution,angular_conus)
import numpy as np
Pi=np.pi
import warnings
warnings.filterwarnings("ignore")
from tkinter import filedialog

M_conf=[{'Mirror_type' : 'spherical', 'p' : 1000, 'q' : 1000, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100],"R": 1000}
        # ,{'Mirror_type' : 'toroid', 'p' : 100, 'q' : 1000, 'alfa' : -85, 'type' : 'focusing', 'size' : [200,100]}
        # ,{'Mirror_type' : 'ellipsoid', 'p' : 2400, 'q' : 700, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
        ]

# M_conf=[{'Mirror_type' : 'toroid', 'p' : 1000, 'q' : 150, 'alfa' : 85, 'type' : 'collimating', 'size' : [200,100]}
#         ,{'Mirror_type' : 'toroid', 'p' : 150, 'q' : 500, 'alfa' : 85, 'type' : 'focusing', 'size' : [200,100]}
#         # ,{'Mirror_type' : 'ellipsoid', 'p' : 2400, 'q' : 700, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
#         ]

# M_conf=[{'Mirror_type' : 'toroid', 'p' : 1000, 'q' : 500, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
#         # ,{'Mirror_type' : 'toroid', 'p' : 150, 'q' : 500, 'alfa' : -85, 'type' : 'focusing', 'size' : [200,100]}
#         # ,{'Mirror_type' : 'ellipsoid', 'p' : 2400, 'q' : 700, 'alfa' : 85, 'type' : 'refocusing', 'size' : [200,100]}
#         ]

#generate rays
# ConAng=1*10**-3 #in rad; angle of the conus
# Nrays=100 # number of rays
# rays=angular_conus(ConAng,Nrays)

AngularWidth=1*10**-3 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
Nrays=1000 # number of rays
rays=gaussian_angular_distribution(AngularWidth,Nrays)

# AngularWidth=1*10**-3 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
# Wradius=0.06 #starting beam radius on the e**-2 intensity level 
# Nrays=4000 # number of rays
# rays=gaussian_spaceYangle_distribution(Wradius,int(np.sqrt(Nrays)),AngularWidth,int(np.sqrt(Nrays)),Nsub0=6)


#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
trace.propagate()
trace.show_rays()
trace.show_plane()