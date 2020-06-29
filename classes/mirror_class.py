"""
mirror class

@author: Slawa

to do:
add check intersection taking the mirror size into account
paraboloid parameters at AOI<45 (wrong not collimating)
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)

from scipy.optimize import fsolve
# import scipy
# from ray_class import ray_class
import sympy
import numpy as np
Pi=np.pi
import error_class as ER
from intersect import cross_spher, cross_plane

class mirror_class():
    def __init__(self,X0,XM0,p=0,q=0,A=[0,0,0],R=0,r=0,f=0,Mirror_type='toroid',angle=0):
        """X0 - centre of the mirror surface (for example sphere centre of a spherical mirror); 
        XM0 centre of the mirror; p, q input and output arm lengths
        A vector defining orientation. only for plane mirror
        R, r - radii of a surface. spherical mirror requires R; toroid and ellipsoid both; R>=r
        angle is required for ellipsoid orientation
        available mirror types: plane, spherical, ellipsoid, toroid"""
        
        if Mirror_type=='spherical':
            self.Type='spherical'
            self.F_surface=spher
            self.n2surface=n_spher
            self.R=R
            self.X0=X0
            self.Args=(X0,R)
            #intersection function
            # x,y,z,xm0,ym0,zm0,amx,amy,amz,xb0,yb0,zb0,abx,aby,abz,R,r,a,b,c = sympy.symbols('x,y,z,xm0,ym0,zm0,amx,amy,amz,xb0,yb0,zb0,abx,aby,abz,R,r,a,b,c')
            # B2SM=sympy.solve([spher(x,y,z,[xm0,ym0,zm0],R),abz*(x-xb0)-abx*(z-zb0),abz*(y-yb0)-aby*(z-zb0),aby*(x-xb0)-abx*(y-yb0)],[x,y,z])
            # self.B2S=B2SM
            self.p=p
            self.q=q
            self.XM0=XM0
        elif Mirror_type=='plane':
            self.Type='plane'
            self.F_surface=plane
            self.X0=X0
            self.A=A
            self.Args=(X0,A)
            self.p=p
            self.q=q
            self.XM0=XM0
        elif Mirror_type=='toroid':
            self.Type='toroid'
            self.F_surface=toroid
            self.n2surface=n_toroid
            self.R=R
            self.r=r
            self.X0=X0
            self.Args=[X0,R,r]
            self.p=p
            self.q=q
            self.XM0=XM0
        elif Mirror_type=='ellipsoid':
            self.Type='ellipsoid'
            self.F_surface=ellipsoid
            self.n2surface=n_ellipsoid
            self.R=R
            self.r=r
            # self.c=c
            self.X0=X0
            self.Args=[X0,R,r,angle]
            self.p=p
            self.q=q
            self.XM0=XM0
            self.angle=angle
        elif Mirror_type=='paraboloid':
            self.Type='paraboloid'
            self.F_surface=paraboloid
            self.n2surface=n_paraboloid
            self.f=f
            self.X0=X0
            self.Args=[X0,f,angle]
            self.p=p
            self.q=q
            self.XM0=XM0
            self.angle=angle
        elif Mirror_type=='cylinder_z':
            self.Type='cylinder_z'
            self.F_surface=cylinder_z
            self.n2surface=n_cylinder_z
            self.R=R
            self.X0=X0
            self.Args=(X0,R)
            self.p=p
            self.q=q
            self.XM0=XM0
        elif Mirror_type=='plane_ellips_z':
            self.Type='plane_ellips_z'
            self.F_surface=plane_ellips_z
            self.n2surface=n_plane_ellips_z
            self.R=R
            self.r=r
            # self.c=c
            self.X0=X0
            self.Args=[X0,R,r,angle]
            self.p=p
            self.q=q
            self.XM0=XM0
            self.angle=angle
    
    Type='toroid'
    Dxy=10 #size
    Dz=10
    R=0 #curvature
    r=0 #second (smaller) curvature for toroid and ellipsoid
    f=0 #focus for paraboloid
    # a=0 #curvature for ellipsoid
    # b=0 #curvature for ellipsoid
    # c=0 #curvature for ellipsoid
    F=0 #focal length
    X0=[0,0,0] #center (of the entire surface)
    A=[0,0,0] #vector defining the tilt of a mirror
    p=0 #input arm length
    q=0 #output arm length
    XM0=[0,0,0] #center of the mirror
    angle=0 #tilt angle for ellipsoid
    F_surface=None #surface function
    n2surface=None #function for the surface normal vector
    error_threshold=10**-6 #threshold to warn that a point is not on a surface
    Args=None
    B2S=None #cross function
    
    def n_vector(self,X1):
        """returns normal vector to the mirror at position X1"""
        if self.Type == 'plane':
            return self.A
        else:
            # print(self.Args)
            return self.n2surface(*X1,*self.Args)
        
    def intersection(self,ray):
        """returns intersection coordinates between a ray and the mirror"""
        rX=ray.X
        rA=ray.A
        Im0th=10**-8 #effective zero im part (to check for the intersection)
        if self.Type == 'spherical':
            # Xcros=cross_spher(self.B2S,self.X0,self.R,rX,rA)
            Xcros=np.array(cross_spher(self.X0,self.R,rX,rA))
            Xs0=Xcros[0]
            # if any([np.abs(sympy.im(x)) > Im0th for x in Xs0]):
            #     #add single point solution check
            #     ray.Lost=True
            #     print('no intersection')
            #     raise ER.SL_exception('no intersection')
            # else:
            #     As0=np.array(Xs0-rX)
            #     normal=self.n_vector(Xs0)
            #     #take the solution in the right direction
            #     if np.sum(rA*As0) > 0 and np.sum(normal*rA) < 0:
            #         return np.array([float(Xcros[0][0]),float(Xcros[0][1]),float(Xcros[0][2])])
            #     else:
            #         return np.array([float(Xcros[1][0]),float(Xcros[1][1]),float(Xcros[1][2])])
            if any([np.imag(x) > Im0th for x in Xs0]):
                #add single point solution check
                ray.Lost=True
                print('no intersection')
                raise ER.SL_exception('no intersection')
            else:
                As0=np.array(Xs0-rX)
                normal=self.n_vector(Xs0)
                #take the solution in the right direction
                if np.sum(rA*As0) > 0 and np.sum(normal*rA) < 0:
                    return Xcros[0]
                else:
                    return Xcros[1]
            
        elif self.Type == 'plane':
            Xcros=cross_plane(self.X0,self.A,rX,rA)
            As0=np.array(Xcros-rX)
            if np.sum(rA*As0) > 0:
                # return np.array([float(Xcros[0].evalf()),float(Xcros[1].evalf()),float(Xcros[2].evalf())])
                return np.array(Xcros)
            else:
                ray.Lost=True
                raise ER.SL_exception('no intersection')
            #add check of real crossing
        
        elif self.Type == 'ellipsoid':
            Xcros, err=cross_ellipsoid(self.X0,self.XM0,self.R,self.r,self.angle,rX,rA)
            errThr=120 #might make sense to increase to ~100
            if np.abs(err) < errThr:
                return np.array(Xcros)
            else:
                ray.Lost=True
                print('no intersection')
                raise ER.SL_exception('no intersecion')
        elif self.Type == 'plane_ellips_z':
            Xcros, err=cross_plane_ellips(self.X0,self.XM0,self.R,self.r,self.angle,rX,rA)
            errThr=120 #might make sense to increase to ~100
            if np.abs(err) < errThr:
                return np.array(Xcros)
            else:
                ray.Lost=True
                print('no intersection')
                raise ER.SL_exception('no intersecion')
            
        elif self.Type == 'toroid':
            Xcros, err=cross_toroid(self.X0,self.XM0,self.R,self.r,rX,rA)
            # print(Xcros, err)
            errThr=np.sum(Xcros**2)**0.5 #might make sense to increase to ~100
            if np.abs(err) < errThr:
                return np.array(Xcros)
            else:
                ray.Lost=True
                print('no intersection')
                raise ER.SL_exception('no intersection')
        elif self.Type == 'paraboloid':
            Xcros, err=cross_paraboloid(self.X0,self.XM0,self.f,self.angle,rX,rA)
            errThr=120
            if np.abs(err) < errThr:
                return np.array(Xcros)
            else:
                ray.Lost=True
                print('no intersection')
                raise ER.SL_exception('no intersecion')
                
        elif self.Type == 'cylinder_z':
            Xcros, err=cross_cylinder_z(self.X0,self.XM0,self.R,rX,rA)
            errThr=120
            if np.abs(err) < errThr:
                return np.array(Xcros)
            else:
                ray.Lost=True
                print('no intersection')
                raise ER.SL_exception('no intersecion')
        
                
        else:
            print('unknown surface')
            raise ER.SL_exception('unknown surface')
            
    
    def reflection(self,ray):
        """reflects a ray from the mirror"""
        if not ray.Lost:
            try:
                X1=self.intersection(ray) #intersection point
            except ER.SL_exception:
                pass
            else:
                n=np.array(self.n_vector(X1))
                n=n/np.sum(n**2)**0.5 #normalize
                a0=ray.A
                #new angle of the ray after reflection (according to Fresnel law)
                a=a0-2*n*np.sum(n*a0) #if reflection from outside of the surface (convex mirror) then sing needs to be opposite
                ray.reflection(X1,a)
        return ray
    
    def remove_complex(self,a):
        """to remove complex numbers from sympy like data"""
        Im0th=10**-15 #effective zero im part
        if sympy.im(a) > Im0th:
            return None
        else:
            return sympy.re(a)
    
    def surf_coordinates(self,X,given = 'xz'):
        """returns coordinates to plot the mirror"""
        #to do: get rid of none and 0
        if given == 'xz':
            """computes y from given x and z coordinates"""
            x,y,z = sympy.symbols('x,y,z')
            Ym_func=sympy.solve(self.F_surface(x,y,z,*self.Args),y)
            Xm=X[0]
            Zm=X[1]
            Ym=[None for i in range(len(Ym_func))]
            indDel=[None for i in range(len(Ym_func))]
            for k in range(len(Ym_func)):
                Ym[k]=[self.remove_complex(Ym_func[k].subs([(x,Xm[i]),(z,Zm[i])]).evalf()) for i in range(len(Xm))]
                indDel[k]=[y == None for y in Ym[k]]
            Out = [[Xm,ym,Zm] for ym in Ym]
            for k in range(len(Out)):
                X=np.delete(Out[k][0],indDel[k])
                Y=np.delete(Out[k][1],indDel[k])
                Z=np.delete(Out[k][2],indDel[k])
                Out[k]=[X,Y,Z]
            return Out
        
        elif given == 'yz':
            """computes x from given y and z coordinates"""
            x,y,z = sympy.symbols('x,y,z')
            Xm_func=sympy.solve(self.F_surface(x,y,z,*self.Args),x)
            Ym=X[0]
            Zm=X[1]
            Xm=[None for i in range(len(Xm_func))]
            indDel=[None for i in range(len(Xm_func))]
            for k in range(len(Xm_func)):
                Xm[k]=[self.remove_complex(Xm_func[k].subs([(y,Ym[i]),(z,Zm[i])]).evalf()) for i in range(len(Ym))]
                indDel[k]=[not x == None for x in Xm[k]]
            Out = [[xm,Ym,Zm] for xm in Xm]
            for k in range(len(Out)):
                X=np.array(Out[k][0])[indDel[k]]
                Y=np.array(Out[k][1])[indDel[k]]
                Z=np.array(Out[k][2])[indDel[k]]
                Out[k]=[X,Y,Z]
            return Out
    
    #================== misalignments=========
    def tiltx(self,phi):
        """for misalignment purpose
        phi tilt in the reflection plane in radians
        the mirror is assumed to be originally in xy plane"""
        if self.Type == 'toroid':
            X0=self.X0 #surface center position
            XM0=self.XM0 #mirror center position
            V=XM0-X0
            V2=np.array([V[0]*np.cos(phi)-V[1]*np.sin(phi),
                         V[0]*np.sin(phi)+V[1]*np.cos(phi),
                         V[2]])
            X02=XM0-V2
            self.X0=X02
            self.Args[0]=X02
            #find new outout arm length
            # R=self.R+self.r
            # f=R*np.cos(angle0+phi)/2
            # p=self.p
            # q2=1/(1/f-1/p)
            # self.q=q2
        elif self.Type == 'ellipsoid' or self.Type == 'plane_ellips_z':
            X0=self.X0 #surface center position
            XM0=self.XM0 #mirror center position
            V=XM0-X0
            V2=np.array([V[0]*np.cos(phi)-V[1]*np.sin(phi),
                         V[0]*np.sin(phi)+V[1]*np.cos(phi),
                         V[2]])
            X02=XM0-V2
            self.X0=X02
            self.Args[0]=X02
            self.angle+=phi
            self.Args[3]=self.angle
        elif self.Type == 'paraboloid':
            X0=self.X0 #surface center position
            XM0=self.XM0 #mirror center position
            V=XM0-X0
            V2=np.array([V[0]*np.cos(phi)-V[1]*np.sin(phi),
                         V[0]*np.sin(phi)+V[1]*np.cos(phi),
                         V[2]])
            X02=XM0-V2
            self.X0=X02
            # print(X0,X02)
            self.Args[0]=X02
            self.angle+=phi
            self.Args[2]=self.angle
        else:
            print("unknown mirror type")
            
    def shift(self,X):
        """shift by xector X"""
        self.X0+=X
        self.XM0+=X


"""equations for surfaces"""
def spher(x,y,z,X0,R):
    """spherical surface equation 
    X0: 3d coordinate of the sphere center, R is the radius"""
    return R**2-(x-X0[0])**2-(y-X0[1])**2-(z-X0[2])**2

def n_spher(x,y,z,X0,R):
    """vector orthogonal to spherical surface"""
    x0,y0,z0=X0
    if R < 0:
        K = -1
    else:
        K=1
    return K*np.array((-2*x + 2*x0, -2*y + 2*y0, -2*z + 2*z0))

def cylinder_z(x,y,z,X0,R):
    """cylindrical mirror with cylinder axis along z"""
    return R**2-(x-X0[0])**2-(y-X0[1])**2

def n_cylinder_z(x,y,z,X0,R):
    """vector orthogonal to cylinder_z surface"""
    x0,y0,z0=X0
    return (-2*x + 2*x0, -2*y + 2*y0, 0)

def cylinder_xy(x,y,z,X0,R,angle):
    """cylindrical mirror with cylinder axis parallel to xy plane"""
    return R**2-((x-X0[0])*np.cos(angle)+(y-X0[1])*np.sin(angle))**2-(z-X0[2])**2

def n_cylinder_xy(x,y,z,X0,R,angle):
    """vector orthogonal to cylinder_xy surface"""
    return (-2*np.cos(angle)*((x-X0[0])*np.cos(angle)+(y-X0[1])*np.sin(angle)), 
            -2*np.sin(angle)*((x-X0[0])*np.cos(angle)+(y-X0[1])*np.sin(angle)),
            -2*(z-X0[2]))

def plane_ellips_z(x,y,z,X0,R,r,angle):
    """cylindrical plane-ellipsoidal mirror with cylinder axis along z"""
    a,b = (R,r)
    return 1-((x-X0[0])*np.cos(angle)+(y-X0[1])*np.sin(angle))**2/a**2-(-(x-X0[0])*np.sin(angle)+(y-X0[1])*np.cos(angle))**2/b**2

def n_plane_ellips_z(x,y,z,X0,R,r,angle):
    """returns vector orthogonal to plane-ellipsoidal mirror in point (x,y,z)"""
    x0,y0,z0=X0
    a,b = (R,r)
    return (-2*np.cos(angle)*((x-x0)*np.cos(angle)+(y-y0)*np.sin(angle))/a**2+2*np.sin(angle)*(-(x-x0)*np.sin(angle)+(y-y0)*np.cos(angle))/b**2, 
            -2*np.sin(angle)*((x-x0)*np.cos(angle)+(y-y0)*np.sin(angle))/a**2-2*np.cos(angle)*(-(x-x0)*np.sin(angle)+(y-y0)*np.cos(angle))/b**2, 
            0)

def ellipsoid(x,y,z,X0,R,r,angle):
    """ellipsoid of revolution
    https://en.wikipedia.org/wiki/Ellipsoid"""
    a,b,c = (R,r,r)
    return 1-((x-X0[0])*np.cos(angle)+(y-X0[1])*np.sin(angle))**2/a**2-(-(x-X0[0])*np.sin(angle)+(y-X0[1])*np.cos(angle))**2/b**2-(z-X0[2])**2/c**2

def n_ellipsoid(x,y,z,X0,R,r,angle):
    """returns vector orthogonal to ellipsoid in point (x,y,z)"""
    x0,y0,z0=X0
    a,b,c = (R,r,r)
    return (-2*np.cos(angle)*((x-x0)*np.cos(angle)+(y-y0)*np.sin(angle))/a**2+2*np.sin(angle)*(-(x-x0)*np.sin(angle)+(y-y0)*np.cos(angle))/b**2, 
            -2*np.sin(angle)*((x-x0)*np.cos(angle)+(y-y0)*np.sin(angle))/a**2-2*np.cos(angle)*(-(x-x0)*np.sin(angle)+(y-y0)*np.cos(angle))/b**2, 
            -2*(z - z0)/c**2)


def toroid(x,y,z,X0,R,r):
    """torus
    R is the distance from the center (X0) of the tube to the center of the torus,
    r is the radius of the tube
    https://en.wikipedia.org/wiki/Torus"""
    return 4*((x-X0[0])**2+(y-X0[1])**2)*R**2-((x-X0[0])**2+(y-X0[1])**2+(z-X0[2])**2+R**2-r**2)**2

def n_toroid(x,y,z,X0,R,r):
    """returns vector orthogonal to toroid in point (x,y,z)"""
    x0,y0,z0=X0
    return -np.array((2*(-R + ((x - x0)**2 + (y - y0)**2)**0.5)*(1.0*x - 1.0*x0)*((x - x0)**2 + (y - y0)**2)**(-0.5), 
            2*(-R + ((x - x0)**2 + (y - y0)**2)**0.5)*(1.0*y - 1.0*y0)*((x - x0)**2 + (y - y0)**2)**(-0.5), 
            2*z - 2*z0))

def plane(x,y,z,X0,N):
    """X0 is the center, N is the normal vector to the plane
    https://en.wikipedia.org/wiki/Plane_(geometry)"""
    return N[0]*(x-X0[0])+N[1]*(y-X0[1])+N[2]*(z-X0[2])

def paraboloid(x,y,z,X0,f,angle):
    """paraboloid of revolution
    f - parent focal length"""
    return ((x-X0[0])*np.cos(angle)+(y-X0[1])*np.sin(angle))-(-(x-X0[0])*np.sin(angle)+(y-X0[1])*np.cos(angle))**2/4/f-(z-X0[2])**2/4/f

def n_paraboloid(x,y,z,X0,f,angle):
    """returns vector orthogonal to toroid in point (x,y,z)"""
    return np.array((np.cos(angle)+2*np.sin(angle)*(-(x-X0[0])*np.sin(angle)+(y-X0[1])*np.cos(angle))/4/f, 
                     np.sin(angle)-2*np.cos(angle)*(-(x-X0[0])*np.sin(angle)+(y-X0[1])*np.cos(angle))/4/f,
                     -2*(z-X0[2])/4/f))

"""equations for surface orientations"""

def toroid_parameters(p,q,alfa):
    """returns R and r for a toroid, p - input arm, q - output arm, 
    alfa is the angle of incidence in degree (from toroid normal: grazing is ~90 degree)"""
    return [2/np.cos(alfa/180*Pi)/(1/p+1/q) , 2*np.cos(alfa/180*Pi)/(1/p+1/q)]

def toroid_position(Xb0,AB0,p,alfa0,R,r):
    """defines the position of the toroid center using the beam parameters 
    Xb0 - starting point and AB directional vector. It needs to be the major/central beam lying in the xy plane
    p - length of the input arm
    concave toroid is assumed"""
    #all mirrors presently are assumed to be in one plain
    alfa=alfa0/180*Pi #degree -> rad
    AB=AB0/np.sum(AB0**2)**0.5 #be sure that it is normalized
    MirCen=Xb0+AB*p #center of the mirror surface
    #vector pointing to the center of the torus
    n=-1*np.array([AB[0]*np.cos(alfa)+AB[1]*np.sin(alfa),AB[1]*np.cos(alfa)-AB[0]*np.sin(alfa),0])
    SufCen=MirCen+n*(R+r) #center of the torus
    return SufCen, MirCen

def plain_position(Xb0,AB0,q):
    """computes plain perpendicuar to the major beam with Xb0 and AB0 parameters (starting point and direction)"""
    AB=AB0/np.sum(AB0**2)**0.5 #be sure that it is normalized
    MX0=Xb0+AB*q #plain position
    Mn=-AB #mirror normal
    return (MX0, Mn)

def spher_position(Xb0,AB0,p,alfa0,R):
    """define position of the spher center using the beam parameters 
    Xb0 and AB0 : starting point and direction
    It needs to be the major/central beam lying in the xy plane
    p - length of the input arm
    concave sphere is assumed"""
    alfa=alfa0/180*Pi #degree -> rad
    AB=AB0/np.sum(AB0**2)**0.5 #to be sure that it is normalized
    MirCen=Xb0+AB*p #center of the mirror surface
    #vector pointing to the center of the sphere
    n=-1*np.array([AB[0]*np.cos(alfa)+AB[1]*np.sin(alfa),AB[1]*np.cos(alfa)-AB[0]*np.sin(alfa),0])
    SufCen=MirCen+n*R #center of the sphere
    return SufCen, MirCen

def ellipsoid_parameters(p,q,alfa):
    """p: input arm, q: output arm, alfa: angle of incidence
    applicable only to the refocusing case 
    (which is, probably, the only one that makes sense for an ellipsoid)"""
    a = (p+q)/2
    b = np.sqrt(p*q*(1+np.cos(2*alfa/180*Pi))/2)
    return a, b

def ellipsoid_position(Xb0,AB0,p,q,alfa,R,r):
    """in principle applicable to any case
    though beam is assumed to be in the xy plane"""
    a,b=(R,r)
    ee=np.sqrt(1-b**2/a**2)
    AB=AB0/np.sum(AB0**2) #normalization
    #position of the mirror center
    XM0=Xb0[0]+AB[0]*p
    YM0=Xb0[1]+AB[1]*p
    ZM0=Xb0[2]+AB[2]*p
    #angle between beam and major ellipsoid axis
    el_angle=np.arcsin(q/(2*a*ee)*np.sin(2*alfa/180*Pi))
    #vector pointing toward ellipsoid centre
    AE=np.array([AB[0]*np.cos(el_angle)-AB[1]*np.sin(el_angle),
                AB[0]*np.sin(el_angle)+AB[1]*np.cos(el_angle),0])
    AE=AE/np.sum(AE**2)**0.5
    #coordinated of the ellipsoid centre
    X0=Xb0[0]+AE[0]*a*ee
    Y0=Xb0[1]+AE[1]*a*ee
    Z0=ZM0
    return np.array([X0,Y0,Z0]), np.array([XM0,YM0,ZM0]), el_angle + np.arctan(AB[1]/AB[0])

def paraboloid_parameters(p,alfa):
    """returns paraboloid parameters for a given configuration, 
    p - efective focal length (input arm for collimation and output arm for focusing), 
    alfa is the angle of incidence in degree (from toroid normal: grazing is ~90 degree)"""
    f=p*np.cos(alfa/180*Pi)**2
    return f
    
def paraboloid_position(Xb0,AB0,p,q,alfa0,f,Type='focusing'):
    """beam is assumed to be in the xy plane"""
    alfa=(90-alfa0)/180*Pi #degree -> rad
    AB=AB0/np.sum(AB0**2) #normalization
    #position of the mirror center
    XM0=Xb0[0]+AB[0]*p
    YM0=Xb0[1]+AB[1]*p
    if Type == 'focusing':
        axis=np.array([AB[0],AB[1]])
        if np.abs(alfa0) >= 45:
            if AB[0] > 0:
                angle = Pi + np.arcsin(AB[1])
            else:
                angle =  np.arcsin(AB[1])
        else:
            if AB[0] > 0:
                angle = Pi - np.arcsin(AB[1])
            else:
                angle =  -np.arcsin(AB[1])
        F=q
        DX=F*np.cos(2*alfa)
        DY=F*np.sin(2*alfa)
        DL=(DX**2+DY**2)**0.5
        DV=DL*np.array([axis[0]*np.cos(2*alfa)-axis[1]*np.sin(2*alfa),axis[0]*np.sin(2*alfa)+axis[1]*np.cos(2*alfa)])
        # print(DX,DY,'\n',DV)
        
        X0,Y0 = [XM0,YM0]+DV+f*axis
    else:
        if np.abs(alfa0) >= 45:
            axis=np.array([AB[0]*np.cos(2*alfa)-AB[1]*np.sin(2*alfa),AB[0]*np.sin(2*alfa)+AB[1]*np.cos(2*alfa)])
            if AB[0] > 0:
                angle = np.arcsin(axis[1])
            else:
                angle = Pi + np.arcsin(axis[1])
        else:
            a0=alfa0/180*Pi
            axis=np.array([AB[0]*np.cos(Pi-2*a0)-AB[1]*np.sin(Pi-2*a0),AB[0]*np.sin(Pi-2*a0)+AB[1]*np.cos(Pi-2*a0)])
            if AB[0] > 0:
                angle = Pi - np.arcsin(axis[1])
            else:
                angle = - np.arcsin(axis[1])
        
        X0,Y0 = [Xb0[0],Xb0[1]]-f*axis
    return np.array([X0,Y0,Xb0[2]]), np.array([XM0,YM0,Xb0[2]]), angle

def cylinder_z_parameters(p,q,alfa0):
    """"""
    alfa=(alfa0)/180*Pi #degree -> rad
    f=1/(1/p+1/q)
    R=2*f/np.cos(alfa)
    return R

def cylinder_z_position(Xb0,AB0,p,q,alfa0,R,Type='refocusing'):
    """"""
    alfa=(alfa0)/180*Pi #degree -> rad
    AB=AB0/np.sum(AB0**2) #normalization
    #position of the mirror center
    XM0=Xb0[0]+AB[0]*p
    YM0=Xb0[1]+AB[1]*p
    ZM0=Xb0[2]
    AE=np.array([AB[0]*np.cos(alfa)+AB[1]*np.sin(alfa),
                -AB[0]*np.sin(alfa)+AB[1]*np.cos(alfa),0])
    #coordinated of the centre
    X0=XM0-AE[0]*R
    Y0=YM0-AE[1]*R
    Z0=ZM0
    
    return np.array([X0,Y0,Z0]), np.array([XM0,YM0,ZM0])
    

"""intersection equations"""
#spher
#sympy symbols
# x,y,z,xm0,ym0,zm0,amx,amy,amz,xb0,yb0,zb0,abx,aby,abz,R,r,a,b,c = sympy.symbols('x,y,z,xm0,ym0,zm0,amx,amy,amz,xb0,yb0,zb0,abx,aby,abz,R,r,a,b,c')

# #the line equation is (x-x0)/ax=(y-y0)/ay=(z-z0)/az where (x0,y0,z0) is the starting point of the line (ax,ay,az) vecter defining its direction
# #B2SM=sympy.solve([spher(x,y,z,[xm0,ym0,zm0],R),abz*(x-xb0)-abx*(z-zb0),abz*(y-yb0)-aby*(z-zb0),aby*(x-xb0)-abx*(y-yb0)],[x,y,z])

# def cross_spher(B2SM,XM0,R0,XB0,AB):
#     """returns intersection coordinated with a beam
#     XM0,R mirror parameters (center and radius), XB0,AB beam parameters, 
#     starting point and directional vector"""
#     Out=np.array([[Xx.subs([(xm0,XM0[0]),(ym0,XM0[1]),(zm0,XM0[2]),(R,R0),
#                             (xb0,XB0[0]),(yb0,XB0[1]),(zb0,XB0[2]),
#                            (abx,AB[0]),(aby,AB[1]),(abz,AB[2])]).evalf() 
#                    for Xx in xyz] for xyz in B2SM])
#     return Out

# B2Pl=sympy.solve([plane(x,y,z,[xm0,ym0,zm0],[amx,amy,amz]),abz*(x-xb0)-abx*(z-zb0),abz*(y-yb0)-aby*(z-zb0),aby*(x-xb0)-abx*(y-yb0)],[x,y,z])

# def cross_plane(XM0,AM0,XB0,AB):
#     """XM0,R mirror parameters (center and radius), XB0,AB beam parameters, 
#     start point and directional vector"""
#     EqOut=[B2Pl[x],B2Pl[y],B2Pl[z]]
#     Out=np.array([Xx.subs([(xm0,XM0[0]),(ym0,XM0[1]),(zm0,XM0[2]),
#                             (amx,AM0[0]),(amy,AM0[1]),(amz,AM0[2]),
#                             (xb0,XB0[0]),(yb0,XB0[1]),(zb0,XB0[2]),
#                            (abx,AB[0]),(aby,AB[1]),(abz,AB[2])]).evalf() for Xx in EqOut])
#     return Out
 
"""solutions for toroid and ellipsoid are numerical"""

def cross_toroid(X0,XM0,R0,r0,XB0,AB):
    """X0,XM0,R,r mirror parameters (center and radius), XB0,AB beam parameters, 
    start point and directional vector"""
    
    x0,y0,z0=X0
    R=R0
    r=r0
    xb0,yb0,zb0=XB0
    abx,aby,abz=AB
    indA=AB.argsort()
    Eff0=10**-3 #efective zero for comparisson
    if not indA[0] < Eff0:
        def eq(var):
            x=var
            z=zb0+abz/abx*(x-xb0)
            y=yb0+aby/abx*(x-xb0)
            return toroid(x,y,z,[x0,y0,z0],R,r)
        
        fs=fsolve(eq,XM0[0],factor=0.1,xtol=10**-15)
        Out=np.array([fs[0],yb0+aby/abx*(fs[0]-xb0),zb0+abz/abx*(fs[0]-xb0)])
    elif not indA[1] < Eff0:
        def eq(var):
            y=var
            z=zb0+abz/aby*(y-yb0)
            x=xb0+abx/aby*(y-yb0)
            return toroid(x,y,z,[x0,y0,z0],R,r)
        
        fs=fsolve(eq,XM0[1],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/aby*(fs[0]-yb0),fs[0],zb0+abz/aby*(fs[0]-yb0)])
    else:
        def eq(var):
            z=var
            y=yb0+aby/abz*(z-zb0)
            x=xb0+abx/abz*(z-zb0)
            return toroid(x,y,z,[x0,y0,z0],R,r)
        
        fs=fsolve(eq,XM0[2],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/abz*(fs[0]-zb0),yb0+aby/abz*(fs[0]-zb0),fs[0]])
        
    return Out, eq(fs[0])


def cross_ellipsoid(X0,XM0,R0,r0,angle,XB0,AB):
    """X0,XM0,R,r are mirror parameters (surface and mirror part), XB0,AB beam parameters, 
    start point and directional vector"""
    x0,y0,z0=X0
    R=R0
    r=r0
    xb0,yb0,zb0=XB0
    abx,aby,abz=AB
    indA=AB.argsort()
    Eff0=10**-3 #efective zero for comparisson
    if not indA[0] < Eff0:
        def eq(var):
            x=var
            z=zb0+abz/abx*(x-xb0)
            y=yb0+aby/abx*(x-xb0)
            return ellipsoid(x,y,z,[x0,y0,z0],R,r,angle)
        
        fs=fsolve(eq,XM0[0],factor=0.1,xtol=10**-15)
        Out=np.array([fs[0],yb0+aby/abx*(fs[0]-xb0),zb0+abz/abx*(fs[0]-xb0)])
    elif not indA[1] < Eff0:
        def eq(var):
            y=var
            z=zb0+abz/aby*(y-yb0)
            x=xb0+abx/aby*(y-yb0)
            return ellipsoid(x,y,z,[x0,y0,z0],R,r,angle)
        
        fs=fsolve(eq,XM0[1],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/aby*(fs[0]-yb0),fs[0],zb0+abz/aby*(fs[0]-yb0)])
    else:
        def eq(var):
            z=var
            y=yb0+aby/abz*(z-zb0)
            x=xb0+abx/abz*(z-zb0)
            return ellipsoid(x,y,z,[x0,y0,z0],R,r,angle)
        
        fs=fsolve(eq,XM0[2],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/abz*(fs[0]-zb0),yb0+aby/abz*(fs[0]-zb0),fs[0]])

    return Out, eq(fs[0])

def cross_paraboloid(X0,XM0,f0,angle,XB0,AB):
    """X0,XM0,f are mirror parameters (center (surface and mirror part) and radii), XB0,AB beam parameters, 
    start point and directional vector"""
    x0,y0,z0=X0
    f=f0
    xb0,yb0,zb0=XB0
    abx,aby,abz=AB
    indA=AB.argsort()
    Eff0=10**-3 #efective zero for comparisson
    if not indA[0] < Eff0:
        def eq(var):
            x=var
            z=zb0+abz/abx*(x-xb0)
            y=yb0+aby/abx*(x-xb0)
            return paraboloid(x,y,z,[x0,y0,z0],f,angle)
        
        fs=fsolve(eq,XM0[0],factor=0.1,xtol=10**-15)
        Out=np.array([fs[0],yb0+aby/abx*(fs[0]-xb0),zb0+abz/abx*(fs[0]-xb0)])
    elif not indA[1] < Eff0:
        def eq(var):
            y=var
            z=zb0+abz/aby*(y-yb0)
            x=xb0+abx/aby*(y-yb0)
            return paraboloid(x,y,z,[x0,y0,z0],f,angle)
        
        fs=fsolve(eq,XM0[1],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/aby*(fs[0]-yb0),fs[0],zb0+abz/aby*(fs[0]-yb0)])
    else:
        def eq(var):
            z=var
            y=yb0+aby/abz*(z-zb0)
            x=xb0+abx/abz*(z-zb0)
            return paraboloid(x,y,z,[x0,y0,z0],f,angle)
        
        fs=fsolve(eq,XM0[2],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/abz*(fs[0]-zb0),yb0+aby/abz*(fs[0]-zb0),fs[0]])

    return Out, eq(fs[0])

def cross_cylinder_z(X0,XM0,R,XB0,AB):
    """X0,XM0,R are mirror parameters (center (surface and mirror part) and radii), XB0,AB beam parameters, 
    start point and directional vector"""
    x0,y0,z0=X0
    xb0,yb0,zb0=XB0
    abx,aby,abz=AB
    indA=AB.argsort()
    Eff0=10**-3 #efective zero for comparisson
    if not indA[0] < Eff0:
        def eq(var):
            x=var
            z=zb0+abz/abx*(x-xb0)
            y=yb0+aby/abx*(x-xb0)
            return cylinder_z(x,y,z,[x0,y0,z0],R)
        
        fs=fsolve(eq,XM0[0],factor=0.1,xtol=10**-15)
        Out=np.array([fs[0],yb0+aby/abx*(fs[0]-xb0),zb0+abz/abx*(fs[0]-xb0)])
    elif not indA[1] < Eff0:
        def eq(var):
            y=var
            z=zb0+abz/aby*(y-yb0)
            x=xb0+abx/aby*(y-yb0)
            return cylinder_z(x,y,z,[x0,y0,z0],R)
        
        fs=fsolve(eq,XM0[1],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/aby*(fs[0]-yb0),fs[0],zb0+abz/aby*(fs[0]-yb0)])
    else:
        def eq(var):
            z=var
            y=yb0+aby/abz*(z-zb0)
            x=xb0+abx/abz*(z-zb0)
            return cylinder_z(x,y,z,[x0,y0,z0],R)
        
        fs=fsolve(eq,XM0[2],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/abz*(fs[0]-zb0),yb0+aby/abz*(fs[0]-zb0),fs[0]])

    return Out, eq(fs[0])

def cross_cylinder_xy(X0,XM0,R,angle,XB0,AB):
    """X0,XM0,R are mirror parameters (center (surface and mirror part) and radii), XB0,AB beam parameters, 
    start point and directional vector"""
    x0,y0,z0=X0
    xb0,yb0,zb0=XB0
    abx,aby,abz=AB
    indA=AB.argsort()
    Eff0=10**-3 #efective zero for comparisson
    if not indA[0] < Eff0:
        def eq(var):
            x=var
            z=zb0+abz/abx*(x-xb0)
            y=yb0+aby/abx*(x-xb0)
            return cylinder_xy(x,y,z,[x0,y0,z0],R,angle)
        
        fs=fsolve(eq,XM0[0],factor=0.1,xtol=10**-15)
        Out=np.array([fs[0],yb0+aby/abx*(fs[0]-xb0),zb0+abz/abx*(fs[0]-xb0)])
    elif not indA[1] < Eff0:
        def eq(var):
            y=var
            z=zb0+abz/aby*(y-yb0)
            x=xb0+abx/aby*(y-yb0)
            return cylinder_xy(x,y,z,[x0,y0,z0],R,angle)
        
        fs=fsolve(eq,XM0[1],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/aby*(fs[0]-yb0),fs[0],zb0+abz/aby*(fs[0]-yb0)])
    else:
        def eq(var):
            z=var
            y=yb0+aby/abz*(z-zb0)
            x=xb0+abx/abz*(z-zb0)
            return cylinder_xy(x,y,z,[x0,y0,z0],R,angle)
        
        fs=fsolve(eq,XM0[2],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/abz*(fs[0]-zb0),yb0+aby/abz*(fs[0]-zb0),fs[0]])

    return Out, eq(fs[0])


def cross_plane_ellips(X0,XM0,R0,r0,angle,XB0,AB):
    """X0,XM0,R,r are mirror parameters (surface and mirror part), XB0,AB beam parameters, 
    start point and directional vector"""
    x0,y0,z0=X0
    R=R0
    r=r0
    xb0,yb0,zb0=XB0
    abx,aby,abz=AB
    indA=AB.argsort()
    Eff0=10**-3 #efective zero for comparisson
    if not indA[0] < Eff0:
        def eq(var):
            x=var
            z=zb0+abz/abx*(x-xb0)
            y=yb0+aby/abx*(x-xb0)
            return plane_ellips_z(x,y,z,[x0,y0,z0],R,r,angle)
        
        fs=fsolve(eq,XM0[0],factor=0.1,xtol=10**-15)
        Out=np.array([fs[0],yb0+aby/abx*(fs[0]-xb0),zb0+abz/abx*(fs[0]-xb0)])
    elif not indA[1] < Eff0:
        def eq(var):
            y=var
            z=zb0+abz/aby*(y-yb0)
            x=xb0+abx/aby*(y-yb0)
            return plane_ellips_z(x,y,z,[x0,y0,z0],R,r,angle)
        
        fs=fsolve(eq,XM0[1],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/aby*(fs[0]-yb0),fs[0],zb0+abz/aby*(fs[0]-yb0)])
    else:
        def eq(var):
            z=var
            y=yb0+aby/abz*(z-zb0)
            x=xb0+abx/abz*(z-zb0)
            return plane_ellips_z(x,y,z,[x0,y0,z0],R,r,angle)
        
        fs=fsolve(eq,XM0[2],factor=0.1,xtol=10**-15)
        Out=np.array([xb0+abx/abz*(fs[0]-zb0),yb0+aby/abz*(fs[0]-zb0),fs[0]])

    return Out, eq(fs[0])