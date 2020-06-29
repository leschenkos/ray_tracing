"""
ray tracing class

@author: Slawa


"""

import os
import sys
Path=os.path.dirname((os.path.abspath(__file__)))
sys.path.append(Path)
from ray_class import ray_class
from mirror_class import mirror_class, toroid_parameters, toroid_position, plain_position, spher_position
from mirror_class import ellipsoid_parameters, ellipsoid_position, paraboloid_parameters, paraboloid_position
from mirror_class import cylinder_z_parameters, cylinder_z_position
import matplotlib.pyplot as plt
import numpy as np
Pi=np.pi
from tkinter import filedialog

class ray_trace():
    def __init__(self,rays,mirrors=None,mirrors_config=None,MrayX=[0,0,0],MrayA=[1,0,0]):
        self.rays=rays
        if mirrors==None:
            self.mirrors=[]
        else:
            self.mirrors=mirrors
        self.Mray=ray_class(MrayX,MrayA,0)
        self.XMray=ray_class([MrayX[0]+0,MrayX[1]+self.Xpos,MrayX[2]+0],
                             [MrayA[0]*np.cos(self.XAng)-MrayA[1]*np.sin(self.XAng),
                              MrayA[1]*np.cos(self.XAng)+MrayA[0]*np.sin(self.XAng),
                              MrayA[2]+0],0)
        if not mirrors_config == None:
            self.mirrors_config=mirrors_config
            self.config()

    rays=None #input data
    Mray=None #major/central ray
    XMray=None #reference ray to orient the output plane 
    mirrors=None
    mirrors_config=None #list to generate mirrors
    Nsteps=0
    Step=0
    XAng=5*10**-4 #angle for reference beam (to define the output plane)
    Xpos=10**-7 #position of the reference beam (to define the output plane)
    
    
    def propagate(self):
        self.Nsteps=len(self.mirrors) #number of ray tracing sections
        self.Step=0
        # print('steps',self.Nsteps)
        
        while self.Step < self.Nsteps:
            # print('step', self.Step)
            for ray in self.rays:
                # print(ray.X)
                self.mirrors[self.Step].reflection(ray)
            self.Step+=1
            
    def config(self,skip_align=False,skipmirror=0):
        ref_ray=self.Mray
        i=0
        if skip_align:
            # for r in self.rays:
                # r.reset()
            self.mirrors.pop(-1) #remove the current end plane and find a new one
            # propagate the reference beams
            # self.Mray=ray_class([0,0,0],[1,0,0],0)
            ref_ray.reset()
            # self.XMray=ray_class([0,self.Xpos,0],[np.cos(self.XAng),np.sin(self.XAng),0],0)
            self.XMray.reset()
            #find new output plane
            mir=self.mirrors
            mconf=self.mirrors_config
            for i in range(len(mconf)):
                m=mir[i]
                if i == skipmirror:
                    m.reflection(ref_ray)
                    q=m.q
                elif m.Type == 'paraboloid':
                    p=m.p
                    q=m.q
                    alfa=mconf[i]['alfa']
                    f=m.f
                    X0, XM0, angle = paraboloid_position(ref_ray.X,ref_ray.A,p,q,alfa,f,mconf[i]['type'])
                    Mir=mirror_class(X0,XM0,p=p,q=q,f=f,Mirror_type='paraboloid',angle=angle)
                    Mir.Dxy , Mir.Dz = mconf[i]['size']
                    # print('parameters: Mirror ', i,' paraboloid',' f: ',f, '\n effective saggital radius:', Rr) 
                    m=Mir
                    mir[i]=Mir
                    m.reflection(ref_ray)
                
                Mx, Mn = plain_position(ref_ray.X,ref_ray.A,q)
                # print(Mx, Mn)
                end=mirror_class(Mx,Mx,A=Mn,Mirror_type='plane')
                end.reflection(ref_ray)
                if ref_ray.Lost:
                    print('ref_ray.Lost')
                ref_ray.A=ref_ray.Astory[-2]
                ref_ray.Astory[-1]=ref_ray.Astory[-2]
                
            # for m in self.mirrors:
            #     m.reflection(self.Mray)
            # if self.Mray.Lost:
            #     print('Mray lost')
                #raise add
            # Mx, Mn = plain_position(self.Mray.X,self.Mray.A,self.mirrors[-1].q)
            # end=mirror_class(Mx,Mx,A=Mn,Mirror_type='plane') #new output plane
            # end.reflection(self.Mray)
            # if self.Mray.Lost:
            #     print('ref_ray.Lost')
            # ref_ray.A=ref_ray.Astory[-2]
            # ref_ray.Astory[-1]=ref_ray.Astory[-2]
        
            end.Dxy=100
            end.Dz=100
            self.mirrors.append(end)
            
            
            self.Mray.reset()
            self.XMray.reset()
            #propagate reference beams
            for m in self.mirrors:
                m.reflection(self.Mray)
                m.reflection(self.XMray)
            if self.Mray.Lost:
                print('Mray lost')
                #raise add
            if self.XMray.Lost:
                #try smaller angle
                # self.XMray=ray_class([0,0,0],[np.cos(self.XAng/10),np.sin(self.XAng/100),0],0)
                self.XMray.reset()
                for m in self.mirrors:
                    m.reflection(self.XMray)
                if self.XMray.Lost:
                    print('XMray lost')
                    #raise add
            
        else:
            for mconf in self.mirrors_config:
                p=mconf['p']
                q=mconf['q']
                alfa=mconf['alfa']
                if mconf['Mirror_type'] == 'toroid':
                    #                print(Rm,rm)
                    if mconf['type'] == 'refocusing':
                        Rm, rm = toroid_parameters(p,q,alfa) #as specified for optical toroid
                        R, r = (Rm-rm, rm) #as specified for a torus (in geometry)
                    elif mconf['type'] == 'collimating':
                        Rm, rm = toroid_parameters(p,10**30,alfa)
                        R, r = (Rm-rm, rm)
                    elif mconf['type'] == 'focusing':
                        Rm, rm = toroid_parameters(10**30,q,alfa)
                        R, r = (Rm-rm, rm)
                    elif mconf['type'] == 'fixed':
                        R, r = (mconf['R']-mconf['r'], mconf['r'])
                    else:
                        print('unknown type')
                        #raise ... add
                    print('parameters: Mirror ', i,' toroid',' R: ',R+r,' r: ', r)   
                    X0, XM0=toroid_position(ref_ray.X,ref_ray.A,p,alfa,R,r)
                    # print(X0, XM0)
                    Mir=mirror_class(X0,XM0,p=p,q=q,R=R,r=r,Mirror_type='toroid')
                    Mir.Dxy , Mir.Dz = mconf['size']
                    # print(2)
                    
                elif mconf['Mirror_type'] == 'ellipsoid':
                    if mconf['type'] == 'refocusing':
                        R, r = ellipsoid_parameters(p,q,alfa)
                        # print(R,r)
                    elif mconf['type'] == 'fixed':
                        R, r = (mconf['R'], mconf['r'])
                    else:
                        print('unknown type')
                        #raise ... add
                            
                    X0, XM0, angle=ellipsoid_position(ref_ray.X,ref_ray.A,p,q,alfa,R,r)
                    # print(X0, XM0, angle)
                    # dx=XM0[0]-X0[0]
                    # dy=XM0[1]-X0[1]
                    e=(1-r**2/R**2)**0.5
                    Y00=p*q*np.sin(2*alfa/180*Pi)/2/R/e
                    X00=R*(1-Y00**2/r**2)**0.5
                    #dx*np.cos(angle)+dy*np.sin(angle), -dx*np.sin(angle)+dy*np.cos(angle)
                    print('parameters: Mirror ', i,' ellipsoid',' R: ',R,' r: ', r, 
                          ' off-axis position xm ', X00) 
    
                    Mir=mirror_class(X0,XM0,p=p,q=q,R=R,r=r,Mirror_type='ellipsoid',angle=angle)
                    Mir.Dxy , Mir.Dz = mconf['size']
                    
                elif mconf['Mirror_type'] == 'plane_ellips_z':
                    if mconf['type'] == 'refocusing':
                        R, r = ellipsoid_parameters(p,q,alfa)
                        # print(R,r)
                    elif mconf['type'] == 'fixed':
                        R, r = (mconf['R'], mconf['r'])
                    else:
                        print('unknown type')
                        #raise ... add
                            
                    X0, XM0, angle=ellipsoid_position(ref_ray.X,ref_ray.A,p,q,alfa,R,r)
                    # print(X0, XM0, angle)
                    # dx=XM0[0]-X0[0]
                    # dy=XM0[1]-X0[1]
                    e=(1-r**2/R**2)**0.5
                    Y00=p*q*np.sin(2*alfa/180*Pi)/2/R/e
                    X00=R*(1-Y00**2/r**2)**0.5
                    #dx*np.cos(angle)+dy*np.sin(angle), -dx*np.sin(angle)+dy*np.cos(angle)
                    print('parameters: Mirror ', i,'plane-ellipsoidal mirror',' R: ',R,' r: ', r, 
                          ' off-axis position xm ', X00) 
    
                    Mir=mirror_class(X0,XM0,p=p,q=q,R=R,r=r,Mirror_type='plane_ellips_z',angle=angle)
                    Mir.Dxy , Mir.Dz = mconf['size']
                    
                elif mconf['Mirror_type'] == 'plane':
                    p=mconf['p']
                    q=mconf['q']
                    alfa=mconf['alfa']
                    XM0, A0 = plain_position(ref_ray.X,ref_ray.A,p)
                    MA= np.array([A0[0]*np.cos(alfa/180*Pi)+A0[1]*np.sin(alfa/180*Pi),
                                 -A0[0]*np.sin(alfa/180*Pi)+A0[1]*np.cos(alfa/180*Pi),
                                  A0[2]]) #vector determining mirror orientation
                    Mir=mirror_class(XM0,XM0,A=MA,p=p,q=q,Mirror_type='plane')
                    Mir.Dxy , Mir.Dz = mconf['size']
                    
                elif mconf['Mirror_type'] == 'spherical':
                    p=mconf['p']
                    q=mconf['q']
                    alfa=mconf['alfa']
                    R=mconf['R']
                    X0, XM0 = spher_position(ref_ray.X,ref_ray.A,p,alfa,R)
                    Mir=mirror_class(X0,XM0,p=p,q=q,R=R,Mirror_type='spherical')
                    Mir.Dxy , Mir.Dz = mconf['size']
                    
                elif mconf['Mirror_type'] == 'paraboloid':
                    p=mconf['p']
                    q=mconf['q']
                    alfa=mconf['alfa']
                    if mconf['type'] == 'collimating':
                        f=paraboloid_parameters(p,alfa)
                        Rr=p*np.abs(np.sin(2*(Pi/2-alfa/180*Pi)))
                    else:
                        f=paraboloid_parameters(q,alfa)
                        Rr=q*np.abs(np.sin(2*(Pi/2-alfa/180*Pi)))
                    X0, XM0, angle = paraboloid_position(ref_ray.X,ref_ray.A,p,q,alfa,f,mconf['type'])
                    Mir=mirror_class(X0,XM0,p=p,q=q,f=f,Mirror_type='paraboloid',angle=angle)
                    Mir.Dxy , Mir.Dz = mconf['size']
                    print('parameters: Mirror ', i,' paraboloid',' f: ',f, '\n effective saggital radius:', Rr) 
                
                elif mconf['Mirror_type'] == 'cylinder_z':
                    p=mconf['p']
                    q=mconf['q']
                    alfa=mconf['alfa']
                    R=cylinder_z_parameters(p,q,alfa) #refocusing
                    X0, XM0 = cylinder_z_position(ref_ray.X,ref_ray.A,p,q,alfa,R)
                    # print(X0, XM0)
                    Mir=mirror_class(X0,XM0,p=p,q=q,R=R,Mirror_type='cylinder_z')
                    Mir.Dxy , Mir.Dz = mconf['size']
                    print('parameters: Mirror ', i,' cylinder_z',' R: ',R)
                
                else:
                    print('unknown mirror type')
                    #raise ... add
                    
                self.mirrors.append(Mir)
                Mir.reflection(ref_ray)
                i+=1
                #output plane
                Mx, Mn = plain_position(ref_ray.X,ref_ray.A,q)
                # print(Mx, Mn)
                end=mirror_class(Mx,Mx,A=Mn,Mirror_type='plane')
                end.reflection(ref_ray)
                if ref_ray.Lost:
                    print('ref_ray.Lost')
                ref_ray.A=ref_ray.Astory[-2]
                ref_ray.Astory[-1]=ref_ray.Astory[-2]
        
            end.Dxy=100
            end.Dz=100
            self.mirrors.append(end) #finish at the out plane of the last mirror
            
            # self.Mray=ray_class([0,0,0],[1,0,0],0)
            self.Mray.reset()
            #propagate reference beams
            for m in self.mirrors:
                m.reflection(self.Mray)
                m.reflection(self.XMray)
            if self.Mray.Lost:
                print('Mray lost')
                #raise add
            if self.XMray.Lost:
                #try smaller angle
                # self.XMray=ray_class([0,0,0],[np.cos(self.XAng/10),np.sin(self.XAng/100),0],0)
                self.XMray.reset()
                for m in self.mirrors:
                    m.reflection(self.XMray)
                if self.XMray.Lost:
                    print('XMray lost')
                    #raise add
    
    def find_focus(self,collimated=False):
        """finds focus (tangential) of the system"""
        if collimated:
            #reference rays 
            X0=self.Mray.Xstory[0]
            A0=self.Mray.Astory[0]
            # print(X0,A0)
            # print('prop')
            r1=ray_class([X0[0],X0[1]+1,X0[2]],A0)
            r2=ray_class(X0,A0) #Mray
        else:
            #reference rays 
            X0=self.Mray.Xstory[0]
            A0=self.Mray.Astory[0]
            # print(X0,A0)
            # print('prop')
            r1=ray_class(X0, [A0[0]*np.cos(self.XAng)-A0[1]*np.sin(self.XAng),
                              A0[1]*np.cos(self.XAng)+A0[0]*np.sin(self.XAng),
                             A0[2]+0])
            r2=ray_class(X0,A0) #Mray
        #propagate the rays
        self.Nsteps=len(self.mirrors)-1 #exclude the output plane
        self.Step=0
        
        while self.Step < self.Nsteps:
            self.mirrors[self.Step].reflection(r1)
            self.mirrors[self.Step].reflection(r2)
            self.Step+=1
        
        # print(r2.Xstory,r1.Xstory)
        # print(r2.Astory,r1.Astory)
        #find intersection (in xy projection)
        X1=r1.X
        A1=r1.A
        a1=A1[1]/A1[0]
        X2=r2.X
        A2=r2.A
        a2=A2[1]/A2[0]
        Xi=(X2[1]-X1[1]-X2[0]*a2+X1[0]*a1)/(a1-a2)
        Yi=X2[1]+a2*(Xi-X2[0])
        Zi=X2[2]+A2[2]/A2[0]*(Xi-X2[0])
        #output arm length
        q=np.sum((X2-np.array([Xi,Yi,Zi]))**2)**0.5
        return q

    def show_rays(self):
        """show rays and mirrors in the xy plane"""
        rays=self.rays
        mirrors=self.mirrors
        #transmitted beams
        indr=[not r.Lost for r in rays]
        Rays=np.array(rays)[indr]
        #limits
        XYextras=0
        Xx=np.array([np.array(r.Xstory)[:,0] for r in Rays])
        # print(Xx)
        Xmin=Xx.min()
        Xmax=Xx.max()
        Yx=np.array([np.array(r.Xstory)[:,1] for r in Rays])
        Ymin=Yx.min()
        Ymax=Yx.max()
        for m in mirrors:
            Xm0=np.linspace(float(Xmin)-XYextras,float(Xmax)+XYextras,500)
            Zm0=np.zeros(len(Xm0))
            Ym0=m.surf_coordinates([Xm0,Zm0],given = 'xz')
            if not Ym0 == []:
                Xm=Ym0[0][0]
                Ym=Ym0[0][1]
                Zm=Ym0[0][1]
                
                for k in range(len(Ym0)-1):
                    Xm=np.concatenate((Xm,Ym0[k+1][0]), axis=0)
                    Zm=np.concatenate((Zm,Ym0[k+1][2]), axis=0)
                    Ym=np.concatenate((Ym,Ym0[k+1][1]), axis=0)
                
                Xmc=[]
                Ymc=[]
                X0,Y0,Z0 = m.XM0
                for k in range(len(Xm)):
                    if not Ym[k] == None:
                        if ((Xm[k]-X0)**2+(Ym[k]-Y0)**2)**0.5 < m.Dxy/2:
                            Xmc.append(Xm[k])
                            Ymc.append(Ym[k])
                
                line,=plt.plot(Xmc,Ymc,linestyle='', marker='o', markersize=0.7)
                line.set_label('mirror')
        
        for r in rays:
            path=np.array(r.Xstory)
            plt.plot(path[:,0],path[:,1])
        plt.axes().set_aspect('equal', 'datalim')
        plt.ylim(float(Ymin)-XYextras,float(Ymax)+XYextras)
        plt.title('Top view')
        plt.xlabel('(mm)')
        plt.ylabel('(mm)')
        plt.legend(tuple(['Mirror' + str(i) for i in range(len(mirrors)-1)] + ['output plane']),
            loc='upper right')
        plt.show()
        
    def get_plane(self,surface=-1):
        """computes the beam profile (points of rays crossing the output plane)"""
        Xr=[]
        Yr=[]
        #vectors to project to
        Ax0=np.array(self.XMray.Xstory[surface]-self.Mray.Xstory[surface])
        Ax=Ax0/np.sum(Ax0**2)**0.5 #normalize
        Ay0=np.cross(self.Mray.Astory[surface],Ax)
        Ay=Ay0/np.sum(Ay0**2)**0.5 #normalize
        for ray1 in self.rays:
                if not ray1.Lost:
                    Xxr=ray1.Xstory[surface]-self.Mray.Xstory[surface]
                    Xr.append(np.dot(Xxr,Ax))
                    Yr.append(np.dot(Xxr,Ay))
        return np.array(Xr),np.array(Yr)
    
    def show_plane(self,surface=-1):
        """shows beam profile (points of rays crossing the output plane)
        in the plane corresponding to surface (-1 : default - output plane;)
        options surface=n will show the profile on the n-th mirror, -2 is the last mirror"""
        Xr, Yr = self.get_plane(surface)

        plt.plot(Xr*10**3,Yr*10**3,linestyle='', marker='o')
        # plt.axes().set_aspect('equal', 'datalim')
        plt.title('Output plane')
        plt.xlabel('x (um)')
        plt.ylabel('y (um)')
        plt.show()
        print('Number of lost rays: ', np.sum([r.Lost for r in self.rays]))
        Xr=np.array(Xr)
        Yr=np.array(Yr)
        rmsX=Xr.std()
        rmsY=Yr.std()
        xmean=Xr.mean()
        ymean=Yr.mean()
        dist=np.array([((Xr[i]-xmean)**2+(Yr[i]-ymean)**2)**0.5 for i in range(len(Xr))])
        rmsD=dist.mean()
        print('std x ', rmsX, '\t', 'std y ', rmsY,'\n', 'std 2D ', rmsD)
        
    def export_plane(self,surface=-1):
        """exports beam profile (points of rays crossing the output plane)"""
        Xr, Yr = self.get_plane(surface)
        file=filedialog.asksaveasfilename()
        if file != '':
            np.savetxt(file+'.dat', np.concatenate((np.array(Xr).reshape((-1,1))*10**3, np.array(Yr).reshape((-1,1))*10**3),axis=1),
                header='x um \t y um',delimiter='\t',comments='')
        
    def rms(self,surface=-1):
        """returns RMS spread of the beam
        in the form [RMS in x dimention, RMS in Y dimention, 2D RMS (RMS 2D radius)]"""
        
        Xr, Yr = self.get_plane(surface)
        rmsX=Xr.std()
        rmsY=Yr.std()
        xmean=Xr.mean()
        ymean=Yr.mean()
        dist=np.array([((Xr[i]-xmean)**2+(Yr[i]-ymean)**2)**0.5 for i in range(len(Xr))])
        rmsD=dist.mean()
        return [rmsX,rmsY,rmsD]
        
    # def points_mirror(self,surface=-1):
    #     """for test purposes
    #     show the point of crossing the mirror"""
    #     Xr=[]
    #     Yr=[]
    #     for ray1 in self.rays:
    #         if not ray1.Lost:
    #             Xr.append(ray1.Xstory[surface][0])
    #             Yr.append(ray1.Xstory[surface][1])
    #     plt.plot(Xr,Yr,linestyle='', marker='o',)
    #     plt.axes().set_aspect('equal', 'datalim')
        
    #     m=self.mirrors[surface]
    #     Zm0=np.zeros(len(Xr))
    #     Ym0=m.surf_coordinates([Xr,Zm0],given = 'xz')
    #     if not Ym0 == []:
    #         Xm=Ym0[0][0]
    #         Ym=Ym0[0][1]
    #         Zm=Ym0[0][1]
            
    #         # print(Xm, Ym)
            
    #         for k in range(len(Ym0)-1):
    #             Xm=np.concatenate((Xm,Ym0[k+1][0]), axis=0)
    #             Zm=np.concatenate((Zm,Ym0[k+1][2]), axis=0)
    #             Ym=np.concatenate((Ym,Ym0[k+1][1]), axis=0)
            
    #         line,=plt.plot(Xm,Ym,linestyle='', marker='o', markersize=0.7)
    #         line.set_label('mirror')
        
    #     plt.ylim(min(Yr),max(Yr))
    #     plt.xlim(min(Xr),max(Xr))
    #     plt.show()
        
        
    def aberrations(self,lam=1030*10**-6):
        """returns aberrations of a mirror
        lam : wavelength in mm (to convert optical path difference to phase)"""
        surface=-1
        #rays coordinates in the output plane
        Xr=[]
        Yr=[]
        #vectors to project to
        Ax0=np.array(self.XMray.Xstory[surface]-self.Mray.Xstory[surface])
        Ax=Ax0/np.sum(Ax0**2)**0.5 #normalize
        Ay0=np.cross(self.Mray.Astory[surface],Ax)
        Ay=Ay0/np.sum(Ay0**2)**0.5 #normalize
        for ray1 in self.rays:
                if not ray1.Lost:
                    Xxr=ray1.Xstory[surface]-self.Mray.Xstory[surface]
                    Xr.append(np.dot(Xxr,Ax))
                    Yr.append(np.dot(Xxr,Ay))
        
        Path=np.array([r.length() for r in self.rays])-self.Mray.length()
        Phase=Path/lam*2*Pi
        return [Phase,Xr,Yr]