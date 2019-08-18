import numpy as np
import matplotlib.pyplot as plt

class waveguide():
    def __init__(self,amplitude,phase,ntop,ncenter,nbottom,ktop,kcenter,kbottom,frequency,length,width):
        self.A= amplitude
        self.phase=phase
        self.nt = ntop
        self.nc = ncenter
        self.nb = nbottom
        self.kt = ktop
        self.kc = kcenter
        self.kb = kbottom
        self.freq= frequency
        self.l = length
        self.w=width

    def generate_contact_points(self,theta):
        points =[(0,self.w/2 )]
        d=0
        i=0
        while d<self.l:
            if i==0:
                points+=[(self.w/(2*np.tan(theta)),self.w)]
                d+=self.w/(2*np.tan(theta))
            elif points[-1][1]==self.w:
                points+=[( points[-1][0]+ self.w/np.tan(theta), 0)]
                d+=self.w/(np.tan(theta))
            else:
                points+=[( points[-1][0]+ self.w/np.tan(theta), self.w)]
                d+=self.w/(np.tan(theta))
            i+=1
        return points
    def isForward(self,n,theta):
        if n.imag*np.cos(theta) >0:
            return n.imag*np.cos(theta)>0
        else:
            return n.real*np.cos(theta)>0
    def brewsters_angle(self):
        return np.arctan(self.nt,self.nc)
    def refracted_angle(self,theta1,n1,n2): #calculates refracted angle using snells law
        if self.isForward(n2,theta1):
            return np.arcsin(n1/n2*np.sin(theta1))
        else:
            return np.pi-np.arcsin(n1/n2*np.sin(theta1))
    def fresnel_coeffs(self,polarization,theta1,n1,n2):
        theta2 = self.refracted_angle(theta1,n1,n2)
        if polarization == 's':
            r = (n1*np.cos(theta1)-n2*np.cos(theta2))/(n1*np.cos(theta1)+n2*np.cos(theta2))
            t = (2*n1*np.cos(theta1))/(n1*np.cos(theta1)+n2*np.cos(theta2))
            return (r,t)
        elif polarization=='p':
            r = (n2*np.cos(theta1)-n1*np.cos(theta2))/(n2*np.cos(theta1)+n1*np.cos(theta2))
            t = (2*n1*np.cos(theta1))/(n2*np.cos(theta1)+n1*np.cos(theta2))
            return (r,t)
        else:
            print("sorry didn't enter a correct polarization, please refer to docs")
            return "garbage"
    def Transmission(self,polarization,theta):
        Estart= self.A
        Emag =Estart
        points =self.generate_contact_points(theta)
        if polarization=='p':
            prevpoint=points[0]
            for i in range(1,len(points)):
                currentpoint = points[i]
                Emag=self.fresnel_coeffs(polarization,theta,complex(self.nc,self.kc),complex(self.nt,self.kt))[0]*Emag
                prevpoint=currentpoint
        if polarization=='s':
            prevpoint=points[0]
            for i in range(1,len(points)):
                currentpoint = points[i]
                Emag=self.fresnel_coeffs(polarization,theta,complex(self.nc,self.kc),complex(self.nt,self.kt))[0]*Emag
                prevpoint=currentpoint
        return abs((Emag/Estart))**2
    def loss(self,polarization,theta):
        Estart= self.A
        Emag =Estart
        points =self.generate_contact_points(theta)
        if polarization=='p':
            prevpoint=points[0]
            for i in range(1,len(points)):
                currentpoint = points[i]
                Emag=self.fresnel_coeffs(polarization,theta,complex(self.nc,self.kc),complex(self.nt,self.kt))[0]*Emag
                prevpoint=currentpoint
        if polarization=='s':
            prevpoint=points[0]
            for i in range(1,len(points)):
                currentpoint = points[i]

                Emag=self.fresnel_coeffs(polarization,theta,complex(self.nc,self.kc),complex(self.nt,self.kt))[0]*Emag

                prevpoint=currentpoint
        return .1*np.log10(abs(Estart/Emag)**2)
