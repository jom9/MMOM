import numpy as np
import matplotlib as plt
class multilayers:
    def __init__(self,indexesofRefractions,imagIndex,thicknesses,lamda0,polarization):
        self.P= polarization
        self.lamdavac = lamda0
        self.n = np.zeros((thicknesses.size+2),dtype='complex')+1
        for i in range(1,thicknesses.size+1):
            self.n[i]+= complex(indexesofRefractions[i-1],imagIndex[i-1])
        self.d=thicknesses #thicknesses of each layer
        self.numofLayers= self.n.size
        self.thetas=np.zeros(self.n.size,dtype='complex')
    def isForward(self,n,theta):
        if n.imag*np.cos(theta) >0:
            return n.imag*np.cos(theta)>0
        else:
            return n.real*np.cos(theta)>0
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
    def generate_Angles(self,theta0):
        self.thetas[0]=theta0
        for i in range(1,self.thetas.size):

            self.thetas[i]=self.refracted_angle(self.thetas[i-1],self.n[i-1],self.n[i])
    def Transfer_Matrix(self,i):
        M = np.zeros((2,2),dtype='complex')
        kz = 2*self.n[i]*pi/self.lamdavac*np.cos(self.n[i])
        delta = kz*self.d[i]
        M[0][0]=complex(0,-1)*delta
        M[1][1]=complex(0,1)*delta
        A= np.eye((2,2))
        fc=self.fresnel_coeffs(self.P, self.thetas[i],self.n[i], self.n[i+1])
        A[0][1]=fc[0]
        A[1][0]=fc[0]
        print(A)
        return np.matmul(M,A)/fc[1]
    def global_Transfer_Matrix(self,theta0):
        self.generate_Angles(theta0)
        G_M=np.identity((2),dtype='complex')
        fc=self.fresnel_coeffs(self.P, self.thetas[0],self.n[0], self.n[1])
        G_M[0][1]=fc[0]
        G_M[1][0]=fc[0]
        G_M=G_M*fc[1]
        for i in range(1,self.n.size-2):
            G_M = np.matmul(G_M,self.Transfer_Matrix(i))
        return G_M
    def Transmittance(self,theta0):
        t=1/(self.global_Transfer_Matrix(theta0)[0][0])
        if self.P=='s':
            return abs(t**2)* (((self.n[-2]*np.cos(self.thetas[-2])).real) / (self.n[1]*np.cos(self.thetas[1])).real)
        elif pol == 'p':
            return abs(t**2) * (((self.n[-2]*conj(np.cos(self.thetas[-2]))).real) /(self.n[1]*conj(np.cos(self.thetas[1]))).real)

    def Reflectance(self,theta0):
        r=self.global_Transfer_Matrix(theta0)[1][0]/self.global_Transfer_Matrix(theta0)[0][0]
        return abs(r**2)
th=np.array([600])
k = np.array([0])
n = np.array([3.42])
lamda0=61

m=multilayers(n,k,th,lamda0,'s')
print(m.Reflectance(0) ,m.Transmittance(0))
