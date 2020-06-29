"""
analytical intersections

@author: Slawa
"""

import sympy
from numpy import sqrt

def spher(x,y,z,X0,R):
    """spherical surface equation 
    X0: 3d coordinate of the sphere center, R is the radius"""
    return R**2-(x-X0[0])**2-(y-X0[1])**2-(z-X0[2])**2

def plane(x,y,z,X0,N):
    """X0 is the center, N is the normal vector to the plane
    https://en.wikipedia.org/wiki/Plane_(geometry)"""
    return N[0]*(x-X0[0])+N[1]*(y-X0[1])+N[2]*(z-X0[2])

def cross_spher(XM0,R,XB0,AB):
    """returns intersection coordinated with a beam
    XM0,R mirror parameters (center and radius), XB0,AB beam parameters, 
    starting point and directional vector"""
    (xm0,ym0,zm0)=XM0
    (xb0,yb0,zb0)=XB0
    (abx,aby,abz)=AB

    return [((abx**2*xm0 - abx*aby*yb0 + abx*aby*ym0 - abx*abz*zb0 + abx*abz*zm0 - abx*sqrt(R**2*abx**2 + R**2*aby**2 + R**2*abz**2 - abx**2*yb0**2 + 2*abx**2*yb0*ym0 - abx**2*ym0**2 - abx**2*zb0**2 + 2*abx**2*zb0*zm0 - abx**2*zm0**2 + 2*abx*aby*xb0*yb0 - 2*abx*aby*xb0*ym0 - 2*abx*aby*xm0*yb0 + 2*abx*aby*xm0*ym0 + 2*abx*abz*xb0*zb0 - 2*abx*abz*xb0*zm0 - 2*abx*abz*xm0*zb0 + 2*abx*abz*xm0*zm0 - aby**2*xb0**2 + 2*aby**2*xb0*xm0 - aby**2*xm0**2 - aby**2*zb0**2 + 2*aby**2*zb0*zm0 - aby**2*zm0**2 + 2*aby*abz*yb0*zb0 - 2*aby*abz*yb0*zm0 - 2*aby*abz*ym0*zb0 + 2*aby*abz*ym0*zm0 - abz**2*xb0**2 + 2*abz**2*xb0*xm0 - abz**2*xm0**2 - abz**2*yb0**2 + 2*abz**2*yb0*ym0 - abz**2*ym0**2) + aby**2*xb0 + abz**2*xb0)/(abx**2 + aby**2 + abz**2),
             (abx**2*yb0 - abx*aby*xb0 + abx*aby*xm0 + aby**2*ym0 - aby*abz*zb0 + aby*abz*zm0 - aby*sqrt(R**2*abx**2 + R**2*aby**2 + R**2*abz**2 - abx**2*yb0**2 + 2*abx**2*yb0*ym0 - abx**2*ym0**2 - abx**2*zb0**2 + 2*abx**2*zb0*zm0 - abx**2*zm0**2 + 2*abx*aby*xb0*yb0 - 2*abx*aby*xb0*ym0 - 2*abx*aby*xm0*yb0 + 2*abx*aby*xm0*ym0 + 2*abx*abz*xb0*zb0 - 2*abx*abz*xb0*zm0 - 2*abx*abz*xm0*zb0 + 2*abx*abz*xm0*zm0 - aby**2*xb0**2 + 2*aby**2*xb0*xm0 - aby**2*xm0**2 - aby**2*zb0**2 + 2*aby**2*zb0*zm0 - aby**2*zm0**2 + 2*aby*abz*yb0*zb0 - 2*aby*abz*yb0*zm0 - 2*aby*abz*ym0*zb0 + 2*aby*abz*ym0*zm0 - abz**2*xb0**2 + 2*abz**2*xb0*xm0 - abz**2*xm0**2 - abz**2*yb0**2 + 2*abz**2*yb0*ym0 - abz**2*ym0**2) + abz**2*yb0)/(abx**2 + aby**2 + abz**2), 
             -abz*sqrt(R**2*abx**2 + R**2*aby**2 + R**2*abz**2 - abx**2*yb0**2 + 2*abx**2*yb0*ym0 - abx**2*ym0**2 - abx**2*zb0**2 + 2*abx**2*zb0*zm0 - abx**2*zm0**2 + 2*abx*aby*xb0*yb0 - 2*abx*aby*xb0*ym0 - 2*abx*aby*xm0*yb0 + 2*abx*aby*xm0*ym0 + 2*abx*abz*xb0*zb0 - 2*abx*abz*xb0*zm0 - 2*abx*abz*xm0*zb0 + 2*abx*abz*xm0*zm0 - aby**2*xb0**2 + 2*aby**2*xb0*xm0 - aby**2*xm0**2 - aby**2*zb0**2 + 2*aby**2*zb0*zm0 - aby**2*zm0**2 + 2*aby*abz*yb0*zb0 - 2*aby*abz*yb0*zm0 - 2*aby*abz*ym0*zb0 + 2*aby*abz*ym0*zm0 - abz**2*xb0**2 + 2*abz**2*xb0*xm0 - abz**2*xm0**2 - abz**2*yb0**2 + 2*abz**2*yb0*ym0 - abz**2*ym0**2)/(abx**2 + aby**2 + abz**2) + (abx**2*zb0 - abx*abz*xb0 + abx*abz*xm0 + aby**2*zb0 - aby*abz*yb0 + aby*abz*ym0 + abz**2*zm0)/(abx**2 + aby**2 + abz**2)),
            ((abx**2*xm0 - abx*aby*yb0 + abx*aby*ym0 - abx*abz*zb0 + abx*abz*zm0 + abx*sqrt(R**2*abx**2 + R**2*aby**2 + R**2*abz**2 - abx**2*yb0**2 + 2*abx**2*yb0*ym0 - abx**2*ym0**2 - abx**2*zb0**2 + 2*abx**2*zb0*zm0 - abx**2*zm0**2 + 2*abx*aby*xb0*yb0 - 2*abx*aby*xb0*ym0 - 2*abx*aby*xm0*yb0 + 2*abx*aby*xm0*ym0 + 2*abx*abz*xb0*zb0 - 2*abx*abz*xb0*zm0 - 2*abx*abz*xm0*zb0 + 2*abx*abz*xm0*zm0 - aby**2*xb0**2 + 2*aby**2*xb0*xm0 - aby**2*xm0**2 - aby**2*zb0**2 + 2*aby**2*zb0*zm0 - aby**2*zm0**2 + 2*aby*abz*yb0*zb0 - 2*aby*abz*yb0*zm0 - 2*aby*abz*ym0*zb0 + 2*aby*abz*ym0*zm0 - abz**2*xb0**2 + 2*abz**2*xb0*xm0 - abz**2*xm0**2 - abz**2*yb0**2 + 2*abz**2*yb0*ym0 - abz**2*ym0**2) + aby**2*xb0 + abz**2*xb0)/(abx**2 + aby**2 + abz**2), 
             (abx**2*yb0 - abx*aby*xb0 + abx*aby*xm0 + aby**2*ym0 - aby*abz*zb0 + aby*abz*zm0 + aby*sqrt(R**2*abx**2 + R**2*aby**2 + R**2*abz**2 - abx**2*yb0**2 + 2*abx**2*yb0*ym0 - abx**2*ym0**2 - abx**2*zb0**2 + 2*abx**2*zb0*zm0 - abx**2*zm0**2 + 2*abx*aby*xb0*yb0 - 2*abx*aby*xb0*ym0 - 2*abx*aby*xm0*yb0 + 2*abx*aby*xm0*ym0 + 2*abx*abz*xb0*zb0 - 2*abx*abz*xb0*zm0 - 2*abx*abz*xm0*zb0 + 2*abx*abz*xm0*zm0 - aby**2*xb0**2 + 2*aby**2*xb0*xm0 - aby**2*xm0**2 - aby**2*zb0**2 + 2*aby**2*zb0*zm0 - aby**2*zm0**2 + 2*aby*abz*yb0*zb0 - 2*aby*abz*yb0*zm0 - 2*aby*abz*ym0*zb0 + 2*aby*abz*ym0*zm0 - abz**2*xb0**2 + 2*abz**2*xb0*xm0 - abz**2*xm0**2 - abz**2*yb0**2 + 2*abz**2*yb0*ym0 - abz**2*ym0**2) + abz**2*yb0)/(abx**2 + aby**2 + abz**2), 
             abz*sqrt(R**2*abx**2 + R**2*aby**2 + R**2*abz**2 - abx**2*yb0**2 + 2*abx**2*yb0*ym0 - abx**2*ym0**2 - abx**2*zb0**2 + 2*abx**2*zb0*zm0 - abx**2*zm0**2 + 2*abx*aby*xb0*yb0 - 2*abx*aby*xb0*ym0 - 2*abx*aby*xm0*yb0 + 2*abx*aby*xm0*ym0 + 2*abx*abz*xb0*zb0 - 2*abx*abz*xb0*zm0 - 2*abx*abz*xm0*zb0 + 2*abx*abz*xm0*zm0 - aby**2*xb0**2 + 2*aby**2*xb0*xm0 - aby**2*xm0**2 - aby**2*zb0**2 + 2*aby**2*zb0*zm0 - aby**2*zm0**2 + 2*aby*abz*yb0*zb0 - 2*aby*abz*yb0*zm0 - 2*aby*abz*ym0*zb0 + 2*aby*abz*ym0*zm0 - abz**2*xb0**2 + 2*abz**2*xb0*xm0 - abz**2*xm0**2 - abz**2*yb0**2 + 2*abz**2*yb0*ym0 - abz**2*ym0**2)/(abx**2 + aby**2 + abz**2) + (abx**2*zb0 - abx*abz*xb0 + abx*abz*xm0 + aby**2*zb0 - aby*abz*yb0 + aby*abz*ym0 + abz**2*zm0)/(abx**2 + aby**2 + abz**2))]

def cross_plane(XM0,AM0,XB0,AB):
    """XM0,R mirror parameters (center and radius), XB0,AB beam parameters, 
    start point and directional vector"""
    (xm0,ym0,zm0)=XM0
    (amx,amy,amz)=AM0
    (xb0,yb0,zb0)=XB0
    (abx,aby,abz)=AB
    return [(abx*amx*xm0 - abx*amy*yb0 + abx*amy*ym0 - abx*amz*zb0 + abx*amz*zm0 + aby*amy*xb0 + abz*amz*xb0)/(abx*amx + aby*amy + abz*amz),
            (abx*amx*yb0 - aby*amx*xb0 + aby*amx*xm0 + aby*amy*ym0 - aby*amz*zb0 + aby*amz*zm0 + abz*amz*yb0)/(abx*amx + aby*amy + abz*amz),
            (abz*(amx*xm0 + amy*ym0 + amz*zm0) + amx*(abx*zb0 - abz*xb0) + amy*(aby*zb0 - abz*yb0))/(abx*amx + aby*amy + abz*amz)]




# x,y,z,xm0,ym0,zm0,amx,amy,amz,xb0,yb0,zb0,abx,aby,abz,R,r,a,b,c = sympy.symbols('x,y,z,xm0,ym0,zm0,amx,amy,amz,xb0,yb0,zb0,abx,aby,abz,R,r,a,b,c')


# B2SM=sympy.solve([spher(x,y,z,[xm0,ym0,zm0],R),abz*(x-xb0)-abx*(z-zb0),abz*(y-yb0)-aby*(z-zb0),aby*(x-xb0)-abx*(y-yb0)],[x,y,z])

# # print(B2SM)

# B2Pl=sympy.solve([plane(x,y,z,[xm0,ym0,zm0],[amx,amy,amz]),abz*(x-xb0)-abx*(z-zb0),abz*(y-yb0)-aby*(z-zb0),aby*(x-xb0)-abx*(y-yb0)],[x,y,z])

# print(B2Pl)