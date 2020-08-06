# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 19:32:28 2020

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

M_conf=[{'Mirror_type' : 'toroid', 'p' : 100, 'q' : 100, 'alfa' : 80, 'type' : 'refocusing', 'size' : [20,100]}
        #,{'Mirror_type' : 'toroid', 'p' : 1000, 'q' : 50, 'alfa' : 80, 'type' : 'focusing', 'size' : [200,100]}
        ]

# M_conf=[{'Mirror_type' : 'plane', 'p' : 500, 'q' : 250, 'alfa' : 10, 'type' : 'refocusing', 'size' : [200,100]}
#         ,{'Mirror_type' : 'plane', 'p' : 250, 'q' : 250, 'alfa' : -20, 'type' : 'focusing', 'size' : [200,100]}
#         ,{'Mirror_type' : 'plane', 'p' : 250, 'q' : 500, 'alfa' : 30, 'type' : 'focusing', 'size' : [200,100]}
#         ]

ConAng=10*10**-3 #in rad; angle of the conus
Nrays=100 # number of rays
rays=angular_conus(ConAng,Nrays)

# rays=[ray_class([0,0,0],[np.cos(phi),np.sin(phi),0],0) 
#                         for phi in np.linspace(-ConAng,ConAng,Nrays)]

trace=ray_trace(rays,mirrors_config=M_conf)
trace.propagate()
trace.show_rays()
trace.show_plane()
trace.show_plane(-2)
# trace.points_mirror(-2)