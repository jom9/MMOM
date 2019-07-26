import numpy as np
import matplotlib.pyplot as plt
from sympy import *
class multilayers:
    def __init__(self,indexesofRefractions,imagIndex,thicknesses,lamda0,polarization):
        self.P= polarization
        self.lamdavac = lamda0
        self.n = np.zeros((thicknesses.size),dtype='complex')
        for i in range(0,thicknesses.size):
            self.n[i]+= complex(indexesofRefractions[i],imagIndex[i])
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
        kz = 2*self.n[i]*np.pi/self.lamdavac*np.cos(self.thetas[i])
        delta = kz*self.d[i]
        M[0][0]=np.exp(complex(0,-1)*delta)
        M[1][1]=np.exp(complex(0,1)*delta)
        A= np.identity((2),dtype='complex')
        fc=self.fresnel_coeffs(self.P, self.thetas[i],self.n[i], self.n[i+1])
        A[0][1]=fc[0]
        A[1][0]=fc[0]
        return np.matmul(M,A)/fc[1]
    def global_Transfer_Matrix(self,theta0):
        self.generate_Angles(theta0)
        G_M=np.identity((2),dtype='complex')
        fc=self.fresnel_coeffs(self.P, self.thetas[0],self.n[0], self.n[1])
        G_M[0][1]=fc[0]
        G_M[1][0]=fc[0]
        G_M=G_M*1/fc[1]
        for i in range(1,self.n.size-1):
            G_M = np.matmul(G_M,self.Transfer_Matrix(i))
        return G_M
    def Transmittance(self,theta0):
        t=1/(self.global_Transfer_Matrix(theta0)[0][0])
        if self.P=='s':

            return abs(t)**2#* (((self.n[-2]*np.cos(self.thetas[-2])).real) / (self.n[1]*np.cos(self.thetas[1])).real)
        elif pol == 'p':
            return abs(t)**2 #* (((self.n[-2]*conj(np.cos(self.thetas[-2]))).real) /(self.n[1]*conj(np.cos(self.thetas[1]))).real)

    def Reflectance(self,theta0):
        r=self.global_Transfer_Matrix(theta0)[1][0]/self.global_Transfer_Matrix(theta0)[0][0]
        return abs(r)**2

lamda = symbols('lamda')
fittingCoeff = [4,0.6916,0.8506,0.2533,0.4911,0.1705,0.1114,0.1384,-0.08671]
nexp = ((fittingCoeff[0]*(lamda**4))+(fittingCoeff[1]*(lamda**3))+(fittingCoeff[2]*(lamda**2))+(fittingCoeff[3]*(lamda**1))+(fittingCoeff[4]*(lamda**0)))/((lamda**4)+(fittingCoeff[5]*(lamda**3))+(fittingCoeff[6]*(lamda**2))+(fittingCoeff[7]*(lamda**1))+(fittingCoeff[8]*(lamda**0)))
T=[]
R=[]
l=[]
narray=[]
f = lambdify(lamda,nexp)
step=.1
start =1
for i in range(0,100):

    lamda0=i*step+start
    l+=[lamda0]
    narray+=[f(lamda0)]
    th=np.array([np.inf,10,np.inf])
    k = np.array([0,0,0])
    n = np.array([1,f(lamda0),1])
    m=multilayers(n,k,th,lamda0,'s')

    R+=[m.Reflectance(0)]
    T+=[m.Transmittance(0)]
    print(i,R[i],T[i])
plt.plot(l,R,T)
plt.show()
