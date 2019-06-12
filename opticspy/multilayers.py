import numpy as np
class multilayers:
    def __init__(self,indexofRefraction,thicknesses):
        self.n=indexofRefraction
        self.d=thicknesses
        self.numofLayers= self.n.size

    def Transmittance(self,lamda):

        k0= 2*np.pi/lamda #wave number in vacuum
        kL= k0
        kR =k0
        # wave number to the left and right of the stack, treated as vacuum by default
        M=np.zeros((2,2),dtype=complex) #global transfer matrix
        L=0 # total thickness
        for i in range(self.numofLayers):
            Mi=np.zeros( (2,2) ,dtype=complex) # local transfer matrix
            kprime = self.n[i]*k0 # wave number in material
            Mi[0][0]=np.cos(kprime*self.d[i])
            Mi[0][1]=np.complex(0,np.sin(kprime*self.d[i])*1/kprime)
            Mi[1][0]= np.complex(0,-kprime*np.sin(kprime*self.d[i]))
            Mi[1][1]=np.cos(kprime*self.d[i])
            if i==0:
                M=Mi
            else:
                M= np.matmul(M,Mi)
            L+=self.d[i]
        t= np.complex(0,2*kL*np.exp(np.complex(0,-kR*L))*1/np.complex( -M[1][0]+ kL*kR*M[0][1], kR*M[0][0]*M[1][1] ) )
        r= np.complex(M[1][0]+kR*kL*M[0][1],(kL*M[1][1]-kR*M[0][0]))/np.complex(-M[1][0]+kR*kL*M[0][1],(-kL*M[1][1]+kR*M[0][0]))
        T= (np.absolute(t)**2) *kR/kL
        R= np.absolute(r)**2
        return (T,R)
