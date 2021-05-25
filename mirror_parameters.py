"""
ray tracing v 1.1 // 21.04.20
@author: Slawa

mirror parapeters for ellipsoids, toroids and paraboloids (paraboloids are not fully finished)
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
"""mirror parameters for an ellipsoid
"""

P=2250 #mm - input arm length
Q=900 #mm - output arm length
Alfa=85 #degree - angle of incidence

# P=14500 #mm - input arm length
# Q=125 #mm - output arm length
# Alfa=82 #degree - angle of incidence

M_conf=[{'Mirror_type' : 'ellipsoid', 'p' : P, 'q' : Q, 'alfa' : Alfa, 
         'type' : 'refocusing', 'size' : [200,100]}]

#generate rays
ConAng=1*10**-3 #in rad; angle of the conus
Nrays=100 # number of rays
rays=angular_conus(ConAng,Nrays)

#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
# trace.propagate()
# trace.show_rays()
# trace.show_plane()

#%%
"""mirror parameters for a toroid
"""

P=2450 #mm - input arm length
Q=750 #mm - output arm length
Alfa=85 #degree - angle of incidence

#change type to collimating or focusing for a monochromator
M_conf=[{'Mirror_type' : 'toroid', 'p' : P, 'q' : Q, 'alfa' : Alfa, 'type' : 'refocusing', 
         'size' : [200,100]} ]

#generate rays
ConAng=1*10**-3 #in rad; angle of the conus
Nrays=100 # number of rays
rays=angular_conus(ConAng,Nrays)

#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
# trace.propagate()
# trace.show_rays()
# trace.show_plane()

#%%
"""mirror parameters for an off-axis paraboloid
"""

P=600 #mm input arm length
Q=600 #mm output arm length
Alfa=85 #degree angle of incidence

#change type to focusing if it is a focusing paraboloid
M_conf=[{'Mirror_type' : 'paraboloid', 'p' : P, 'q' : Q, 'alfa' : Alfa, 'type' : 'collimating', 
         'size' : [200,100]} ]

#generate rays
ConAng=1*10**-3 #in rad; angle of the conus
Nrays=100 # number of rays
rays=angular_conus(ConAng,Nrays)

#propagate them
trace=ray_trace(rays,mirrors_config=M_conf)
# trace.propagate()
# trace.show_rays()
# trace.show_plane()