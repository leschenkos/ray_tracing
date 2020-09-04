"""
field propagation class

dimentions: mm

@author: Slawa

to do:
    add log plots
"""
import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)
import numpy as np
Pi=np.pi
import matplotlib.pyplot as plt
from Fresnel_Kirchhoff_diffraction_parallel import FK_diffraction, F_2D
from multiprocessing import Pool, cpu_count

from spherical_mirror_aberrations import mirror_aberration, telescope_aberration
from color_maps import plt_cmap

class field():
    def __init__(self,F=None,X=None,lam=1030*10**-6):
        """F complex filed as a 2d numpy array
        lam wavelength in mm
        X: coordinate vector (a square image is assumed)"""
        self.F=F
        self.lam=lam
        self.k=2*Pi/lam
        self.X=X
        
    lam=1030*10**-6
    F=None #field
    Fin=None
    Xin=None
    X=None
    
    def GaussianBeam(self,X,W,R=10**15,lam=0,W_hole=0):
        """field of a gaussian beam
        X: coordinate vector; W: radius on e**-2 level; R: radius of curvature
        lam: wavelength; Whole: central hole diameter"""
        self.X=X
        Nx=len(X)
        XinM=X[:,None]*np.ones(Nx)
        YinM=np.ones(Nx)[:,None]*X
        if lam==0:
            if not self.lam==0:
                lam=self.lam
            else:
                print('lam=0')
                #add raise
        else:
            self.lam=lam
            self.k=2*Pi/lam
        self.F=E_Gaus(XinM,YinM,W,R,lam,W_hole)
        
    def FK_propagation(self,L,p):
        """Fresnel-Kirchhoff diffraction propagation
        L: distance in mm; p: Pool for parallel computations"""
        #prepare
        DXin=self.X[-1]-self.X[0]
        dxin=self.X[1]-self.X[0]
        Xin=self.X
        k=self.k
        DXout=DXin #assuming preserving the image size
        dxout=1/DXin/k*L*2*Pi
        if dxout > dxin:
            dxout = dxin
        Nout=int(DXout/dxout+1)
        # print(Nout)
        Xout=np.linspace(-DXout/2,DXout/2,Nout) #output dimentions
        #propagate
        Dif=FK_diffraction(self.F,Xin,Xin,Xout,Xout,L,p)
        Dif.propagate()
        #assign the output data
        self.Fin=np.copy(self.F)
        self.F=Dif.Eout
        self.Xin=Xin
        self.X=Xout
        
    def add_lens(self,f):
        """adds phase of a thin lens with focus f (in mm)"""
        X=self.X
        Nx=len(X)
        XinM=X[:,None]*np.ones(Nx)
        YinM=np.ones(Nx)[:,None]*X
        k=self.k
        self.F*=np.exp(-1j*k/2/f*(XinM**2+YinM**2))
        
    def add_tilt_mirror_aberration(self,tilt,R):
        """adds aberrations of a tilted mirror; 
        tilt is the tilt angle in degree; R is the radius of curvature in mm (positive is a focusing mirror)"""
        X=self.X
        Ph=mirror_aberration(R,tilt,X,X,self.lam)
        
        #plot
        # plt.pcolormesh(X,X,Ph,cmap=plt_cmap())
        # plt.xticks(fontsize=16)
        # plt.yticks(fontsize=16)
        # plt.xlabel('x (mm)',fontsize=18)
        # plt.ylabel('y (mm)',fontsize=18)
        # plt.title('phase from mirror tilt')
        # plt.show()
        
        self.F*=np.exp(1j*Ph)
        
    def add_telescope_aberration(self,R1,tilt1,R2,tilt2):
        """adds aberrations of a tilted spherical telescope"""
        X=self.X
        Ph=telescope_aberration(R1,tilt1,R2,tilt2,X,X,self.lam)
        self.F*=np.exp(1j*Ph)
        
    def En_in_circle(self,R,X=[0,0]):
        """returns the part of energy in circle with radius R (<=), centered at X"""
        F=self.F
        En=np.sum(np.abs(F)**2) #full energy
        Xn=self.X
        Nx=len(Xn)
        #find points inside the circle
        XinM=Xn[:,None]*np.ones(Nx)
        YinM=np.ones(Nx)[:,None]*Xn
        ind=(XinM-X[0])**2+(YinM-X[1])**2 <= R**2
        
        Enc=np.sum(np.abs(F[ind])**2) #energy in the circle
        return Enc/En
    
    def Fourier(self,f,DXout,N,p):
        """fourier transform, f - focal length, DXout - size of the focus image
        N - array size (number of points along x and y)
        p: Pool for parallel computations"""
        X=self.X
        kx=np.linspace(-DXout/2,DXout/2,N)
        self.Fin=self.FK_propagation
        self.Xin=X
        self.F=F_2D(self.F,X,kx,self.k,f,p)
        self.X=kx
        
    def Inv_Fourier(self,f,p):
        """inverse fourier transform"""
        X=self.X
        Kx=self.Xin
        F=F_2D(self.F,X,-Kx,self.k,f,p)
        self.Fin=self.Fin
        self.Xin=self.X
        self.F=F
        self.X=Kx
    
    def aperture(self,R,Xa=[0,0]):
        """propagate beam through a hard aperture of radius R centered at Xa"""
        X=self.X
        Nx=len(X)
        XinM=X[:,None]*np.ones(Nx)
        YinM=np.ones(Nx)[:,None]*X
        ind=(XinM-Xa[0])**2+(YinM-Xa[1])**2 > R**2
        self.F[ind]=0
    
    def show_I(self,title=None):
        """plots the intensity distribution"""
        plt.rcParams['figure.dpi']= 300
        plt.rcParams['figure.figsize'] = (5, 5)
        plt.pcolormesh(self.X,self.X,np.abs(self.F)**2,cmap=plt_cmap(),shading='nearest')
        # plt.axis('equal')
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('x (mm)',fontsize=18)
        plt.ylabel('y (mm)',fontsize=18)
        if not title == None:
            plt.title(title)
        plt.show()
        
    def peak_fluence(self):
        """returns peak fluence normalized to unity pulse energy
        multiply by energy and get fluence in J/cm**2"""
        dx=self.X[1]-self.X[2]
        F=self.F
        I=np.abs(F)**2
        En=np.sum(I) #full energy
        return np.max(I)/En/dx**2*100
        
    def diameter(self,method='4sigma'):
        """returns beam diameter
        methods: 4sigma; FWHM; e**-2
        returns [diameter X, diameter Y]"""
        F=self.F
        X=self.X
        I=np.abs(F)**2
        Ix=np.sum(I,axis=1)
        Ix=Ix/Ix.max()
        Iy=np.sum(I,axis=0)
        Iy=Iy/Iy.max()
        if method == 'e**-2':
            ind=Ix >= np.exp(-2)
            indx=np.where(ind == True)
            Dx=X[indx[0][-1]]-X[indx[0][0]]
            ind=Iy >= np.exp(-2)
            indy=np.where(ind == True)
            Dy=X[indy[0][-1]]-X[indy[0][0]]
        elif method == 'FWHM':
            ind=Ix >= 0.5
            indx=np.where(ind == True)
            Dx=X[indx[0][-1]]-X[indx[0][0]]
            ind=Iy >= 0.5
            indy=np.where(ind == True)
            Dy=X[indy[0][-1]]-X[indy[0][0]]
        elif method == '4sigma':
            sigmaX=(np.sum(X**2*Ix)/np.sum(Ix))**0.5
            Dx=4*sigmaX
            sigmaY=(np.sum(X**2*Iy)/np.sum(Iy))**0.5
            Dy=4*sigmaY
        else:
            print('unknown beam diameter method')
        return (Dx,Dy)
    
    def phase_rms(self):
        """returns phase distortions (waited with field amplitude profile)
        returns (RMS, PV)"""
        ph=np.angle(self.F)
        I=np.abs(self.F)
        I*=1/I.max()
        phW=ph*I #weighted phase
        PV=phW.max()
        RMS=phW.std()
        return (RMS,PV)
    
    def save(self):
        pass
    
    def load(self):
        pass
        

def E_Gaus(x,y,W,R,lam,Whole=0):
    """field of a gaussian beam"""
    r=(x**2+y**2)**0.5
    k=2*Pi/lam
    if Whole==0:
        return np.exp(-r**2/W**2-1j*k*r**2/2/R)
    else:
        E=np.exp(-r**2/W**2-1j*k*r**2/2/R)
        ind= r <= Whole
        E[ind]=0
        return E
    
# #============test

# Scaling=2**0.5
# W=1.987*Scaling*1.5 #radius e-2
# W_hole=W*0.2
# print('e-2 radius', W)
# lam=1030*10**-6 #wavelength
# L=1500 #distance
# DXin=W*3/3*2
# Nxin=int(46*Scaling**2*3000/L)+1
# print(Nxin)
# Xin=np.linspace(-DXin,DXin,Nxin)

# Fl=field()
# Fl.GausianBeam(Xin,W,lam=lam,W_hole=W_hole)
# Fl.show_I()
# print(Fl.peak_fluence())

# #if cycle is a wudu magic for parallel computations to work in windows
# if __name__ == '__main__':
#     #initiate parallel computing
#     Ncpu=cpu_count()
#     if Ncpu > 1:
#         Npr=Ncpu-1 #number of parallel processes in the Pool
#     else:
#         Npr=1
#     p=Pool(Npr)
#     # print(p._processes)
#     # Fl.FK_propagation(L,p)
#     # Fl.show_I()
#     # print(Fl.En_in_circle(0.9*W_hole))
#     # print(Fl.peak_fluence())
    
#     # tilt=10
#     # f=750
    
#     # Fl.add_lens(f)
#     # Fl.add_tilt_mirror_aberration(tilt,2*f)
    
#     # Fl.FK_propagation(L,p)
#     # Fl.show_I()
#     # print(Fl.En_in_circle(0.98*W_hole))
    
    
#     f=750
#     Fl.Fourier(f,1200*10**-3,Nxin,p)
#     Fl.show_I()
#     # Fl.aperture(100*10**-3)
#     Fl.Inv_Fourier(f,p)
#     Fl.show_I()
#     print(Fl.En_in_circle(0.8*W_hole))
    
#     p.close()
#     p.join()

#========== test telescope aberrations
# Scaling=2**0.5
# W=20/2 #radius e-2
# W_hole=0
# lam=1030*10**-6 #wavelength
# Nxin=101
# DXin=W*3/3*2
# Xin=np.linspace(-DXin,DXin,Nxin)

# Fl=field()
# Fl.GausianBeam(Xin,W,lam=lam,W_hole=W_hole)
# Fl.show_I()
# print(Fl.phase_rms())
# Fl.add_telescope_aberration(1800,4,-600,-4)
# print(Fl.phase_rms())