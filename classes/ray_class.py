"""
ray class
@author: Slawa
"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)

import numpy as np
import matplotlib.pyplot as plt
Pi=np.pi

class ray_class():
    """starting point, directional vector"""
    def __init__(self,X,A,L=0):
        """X0 is the starting point
        A is the vector defining the ray
        L is the optical path"""
        self.L=L
        self.X=np.array(X)
        self.A=np.array(A)/np.sum(np.array(A)**2) #normalization
        self.Xstory=[np.array(X)]
        self.Astory=[np.array(A)/np.sum(np.array(A)**2)]
    
    X=np.array([0.,0.,0.]) #start coordinate
    A=np.array([0.,0.,0.]) #vector defining the ray
    a=np.array([0.,0.,0.]) #direction
    Lost=False #status if missed a mirror
    Xstory=None #the propagation story
    Astory=None
    
    def length(self):
        if len(self.Xstory) == 1:
            return 0
        else:
            return self.L+np.sum([(np.sum([(self.Xstory[i+1]-self.Xstory[i])**2]))**0.5 
                           for i in range(len(self.Xstory)-1)])
    
    def reflection(self,X,A):
        """define new angle A after reflection at X"""
        X1=X
        self.Xstory.append(np.array(X1))
        self.Astory.append(np.array(A))
        self.X=np.array(X1)
        self.A=np.array(A)/np.sum(np.array(A)**2) #normalization
        
    def reset(self):
        self.X=self.Xstory[0]
        self.A=self.Astory[0]
        self.Xstory=[np.copy(self.X)]
        self.Astory=[np.copy(self.A)]
        self.L=0
 
def angular_conus(angle,Nrays):
    """angle is the conus angle in radians; Nrays in the number of rays"""
    rays=[ray_class([0,0,0],[np.cos(angle),np.sin(angle)*np.cos(phi),np.sin(angle)*np.sin(phi)],0) 
                            for phi in np.linspace(-Pi,Pi,Nrays)]
    return rays

def gaussian_angular_distribution(a,N,X0=[0,0,0],YTilt=0,Nsub0=10):
    """a is the angle defining the divergence on the e**-2 intensity level, N number of rays
    Nsub is the angular discretization"""
    Nsub=Nsub0+2
    Levels=np.linspace(1/Nsub/2,1-1/Nsub/2,Nsub-2)
    angles=(-np.log(Levels)*a**2/2)**0.5
    n=int(N/Nsub) #number of rays in each conus
    rays=[]
    for ConAng in angles:
        
        rays+=[ray_class(X0,ytilt([np.cos(ConAng),np.sin(ConAng)*np.cos(phi),np.sin(ConAng)*np.sin(phi)],YTilt),0) 
                        for phi in np.linspace(-Pi,Pi,n)]
        
        # rays+=[ray_class(X0,rotate(A0,ConAng,phi),0) 
        #                 for phi in np.linspace(-Pi,Pi,n)]
    return rays

def gaussian_spaceYangle_distribution(W,Ns,Wa,Na,Nsub0=10):
    """a angle defining the divergence on the e**-2 intensity level, N number of rays
    Nsub is the angular discretization"""
    Nsub=Nsub0+2
    Levels=np.linspace(1/Nsub/2,1-1/Nsub/2,Nsub-2)
    angles=(-np.log(Levels)*Wa**2/2)**0.5
    radii=(-np.log(Levels)*W**2/2)**0.5
    # print(radii)
    ns=int(Ns/Nsub) #number of rays in each conus
    na=int(Na/Nsub) #number of rays in each conus
    rays=[]
    #central point
    for ConAng in angles:
            rays+=[ray_class([0,0,0],[np.cos(ConAng),np.sin(ConAng)*np.cos(phi),np.sin(ConAng)*np.sin(phi)],0) 
                        for phi in np.linspace(-Pi,Pi,na)]
    #all other points
    for r in radii:
        points=[[0,r*np.cos(a),r*np.sin(a)] for a in np.linspace(-Pi,Pi,ns)]
        # print(points,'\n')
        for p in points:
            for ConAng in angles:
                rays+=[ray_class(p,[np.cos(ConAng),np.sin(ConAng)*np.cos(phi),np.sin(ConAng)*np.sin(phi)],0) 
                            for phi in np.linspace(-Pi,Pi,na)]
    return rays

def gaussian_space_distribution(W,N,X0=[0,0,0],A=[1,0,0],Nsub0=10):
    """W beam radius on the e**-2 intensity level, N number of rays
    Nsub is the angular discretization, A direction vector"""
    Nsub=Nsub0+2
    Levels=np.linspace(1/Nsub/2,1-1/Nsub/2,Nsub-2)
    radii=(-np.log(Levels)*W**2/2)**0.5
    # print(radii)
    ns=int(N/Nsub) #number of rays in each conus
    rays=[]
    #central point
    rays+=[ray_class(X0,A)]
    #all other points
    for r in radii:
        points=[[X0[0]+0,X0[1]+r*np.cos(a),X0[2]+r*np.sin(a)] for a in np.linspace(-Pi,Pi,ns)]
        rays+=[ray_class(p,A) for p in points]

    return rays

def beam2aberrations(a,N,X=[1],Y=[1],Btype='rectengular'):
    """round collimared beam for aberration calculation
    a : radius
    N : namber of rays"""
    rays=[]
    # rays+=[ray_class([0,0,0],[1,0,0],0)]
    if len(X)>1 and len(Y)>1:
        for x in X:
                for y in Y:
                    rays+=[ray_class([0,x,y],[1,0,0],0)]
        return rays
    else:
        if Btype == 'round':
            Fil=0.78 #~fill factor
            S=Pi*a**2
            R=(S/N/Pi*Fil)**0.5
            D=2*R #distance bitween rays
            Nr=int(np.ceil(a/D)) #number of radial steps
            for i in range(1,Nr):
                Ri=D*i
                Ni=int(np.ceil(2*Pi*Ri/D))
                for j in range(Ni):
                    Y=Ri*np.cos(2*Pi/Ni*j)
                    Z=Ri*np.sin(2*Pi/Ni*j)
                    rays+=[ray_class([0,Y,Z],[1,0,0],0)]
            return rays
        else:
            Nx=int(np.ceil(N**0.5))
            X=np.linspace(-a,a,Nx)
            Y=np.linspace(-a,a,Nx)
            for i in range(Nx):
                for j in range(Nx):
                    rays+=[ray_class([0,X[i],Y[j]],[1,0,0],0)]
            return rays,X,Y
    

def rotate(A,theta,phi):
    """3d rotation of vector A; theta around z axis; phi around x axis (to fit the simulation geometry)"""
    R1=np.array([[np.cos(theta),-np.sin(theta),0],
                 [np.sin(theta),np.cos(theta),0],
                 [0,0,1]])
    R2=np.array([[1,0,0],
                 [0,np.cos(phi),-np.sin(phi)],
                 [0,np.sin(phi),np.cos(phi)]])
    return np.dot(R2,np.dot(R1,np.array(A).reshape((3,1)))).reshape(3)

def ytilt(A,phi):
    """3d rotation of vector A;phi around y axis """
    R=np.array([[np.cos(phi),0,-np.sin(phi)],
                 [0,1,0],
                 [np.sin(phi),0,np.cos(phi)]])
    return np.dot(R,np.array(A).reshape((3,1))).reshape(3)

# print(rotate([np.cos(1),0,np.sin(1)],0.1,0.1))

# rays0=beam2aberrations(10,1000)
# Y0=[r.X[1] for r in rays0]
# Z0=[r.X[2] for r in rays0]
# plt.plot(Y0,Z0,linestyle='', marker='o')
# plt.show()
    