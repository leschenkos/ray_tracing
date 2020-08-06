# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 17:46:53 2020

@author: VL
"""

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

# M_conf=[{'Mirror_type' : 'toroid', 'p' : 1500, 'q' : 150, 'alfa' : 80, 'type' : 'collimating', 'size' : [200,100], 'R' : 17300, 'r' : 521}
#         ,{'Mirror_type' : 'toroid', 'p' : 150, 'q' : 150, 'alfa' : -80, 'type' : 'focusing', 'size' : [200,100], 'R' : 1730, 'r' : 52.1}
#         ,{'Mirror_type' : 'toroid', 'p' : 660, 'q' : 600, 'alfa' : 80, 'type' : 'refocusing', 'size' : [200,100], 'R' : 3620, 'r' : 109.2}
#         ]

# M_conf=[{'Mirror_type' : 'toroid', 'p' : 1500, 'q' : 150, 'alfa' : 80, 'type' : 'fixed', 'size' : [200,100], 'R' : 17300, 'r' : 521}
#         ,{'Mirror_type' : 'toroid', 'p' : 150, 'q' : 150, 'alfa' : -80, 'type' : 'fixed', 'size' : [200,100], 'R' : 1730, 'r' : 52.1}
#         ,{'Mirror_type' : 'toroid', 'p' : 660, 'q' : 600, 'alfa' : 80, 'type' : 'fixed', 'size' : [200,100], 'R' : 3620, 'r' : 109.2}
#         ]

# M_conf=[{'Mirror_type' : 'toroid', 'p' : 6000, 'q' : 15, 'alfa' : 75, 'type' : 'fixed', 'size' : [200,100], 'R' : 2050, 'r' : 137.2}
#         ,{'Mirror_type' : 'toroid', 'p' : 15, 'q' : 170, 'alfa' : 75, 'type' : 'fixed', 'size' : [200,100], 'R' : 4213, 'r' : 281.8}
#         ]

# M_conf=[{'Mirror_type' : 'toroid', 'p' : 1000, 'q' : 15, 'alfa' : 75, 'type' : 'collimating', 'size' : [200,100], 'R' : 2050, 'r' : 137.2}
#         ,{'Mirror_type' : 'toroid', 'p' : 15, 'q' : 170, 'alfa' : 75, 'type' : 'focusing', 'size' : [200,100], 'R' : 4213, 'r' : 281.8}
        # ]
        
M_conf=[{'Mirror_type' : 'ellipsoid', 'p' : 1500, 'q' : 500, 'alfa' : 80, 'type' : 'refocusing', 'size' : [200,100], 'R' : 17300, 'r' : 521}
        ]
        

#generate rays
# AngularWidth=0.5*10**-3 #in rad; angle (radius) of the gaussian divergence on the e**-2 intensity level 
# Nrays=1000 # number of rays
# rays=gaussian_angular_distribution(AngularWidth,Nrays)

AngularWidth=1.5*10**-3 #in rad; angle (radus) of the gaussian divergence on the e**-2 intensity level 
Wradius=0.02 #starting beam radius on the e**-2 intensity level 
Nrays=4000 # number of rays
rays=gaussian_spaceYangle_distribution(Wradius,int(np.sqrt(Nrays)*2),AngularWidth,int(np.sqrt(Nrays)/2),Nsub0=6)

# ConAng=1*10**-3 #in rad; angle of the conus
# Nrays=100 # number of rays
# rays=angular_conus(ConAng,Nrays)

#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
trace.propagate()
trace.show_rays()
trace.show_plane()