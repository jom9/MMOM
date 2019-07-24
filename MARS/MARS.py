from sympy import *
import numpy as np
import time
from queue import Queue
import threading
import sys
from scipy import constants
numofthreads=1
threadQ=Queue()
class MARS():
    def __init__(self,thickness,diameter,armLength,startphi, angularvelocity,magnetization,step,resolution ,xlim, ylim,zlim, totaltime,timestep):
        self.thickness=thickness
        self.diameter=diameter
        self.armLength=armLength
        self.phi=startphi
        self.startphi=startphi
        self.angularvelocity=angularvelocity
        self.magnetization=magnetization
        self.xlim=xlim
        self.ylim=ylim
        self.zlim=zlim
        self.totaltime=totaltime
        self.timestep=timestep
        self.resolution=resolution
        self.step = step
    def magnitude(self,vector): #function to calculate the magnitude of a vector
        m=0
        for i in range(vector.size):
            m+=vector[i]**2
        return m**.5
    def curl(self,A,symx,symy,symz): #calculates the curl symbolically in cartesian coordinates
        B=[]
        B+=[diff(A[2],symy)-diff(A[1],symz)]
        B+=[-1*(diff(A[2],symx)-diff(A[0],symz))]
        B+=[diff(A[1],symx)-diff(A[0],symy)]

        return B
    def Magnetic_Field_Slice(self,slice_diameter,zprime):
    # thickness is the thickness of the magnet
    # diameter of the magnet
    #armLength is the distance the magnet is away from the center
    #phi is the current angular position
    #magnetization is the magnatic dipole moments per unit area, this is assumed to be coaxial
    #returns the 3d array of the magentic field

        M = self.magnetization*np.array([np.cos(self.phi),np.sin(self.phi),0]) # magnetization vector

        x,y,z,xprime,yprime= symbols('x y z xprime yprime')
        delR = np.array([x-xprime,y-yprime,z-zprime]) # this is the distance vector, with the primed values denoting source points and the unprimed values denoting field values

        integrand =np.cross(M,delR)/(self.magnitude(delR)**3)

        cornerPos= self.sortbyx(self.currentPos(self.thickness,slice_diameter,self.armLength,self.phi)) #gets the four corners needed to take each rectangluar slice

        posFunction = self.postionFunctions(self.thickness,slice_diameter,self.armLength,self.phi) #creates the lines that trace out borders to slice

        AijkX=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        AijkY=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        AijkZ=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))

        if posFunction[0][0]== np.Infinity:
            Ax =0
        elif(posFunction[0][0]<=0):
            Ax = self.numeric_double_Integral(integrand[0],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],self.step,xprime,yprime)
        else:
            Ax = self.numeric_double_Integral(integrand[0],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[0][0]+posFunction[0][1],self.step,xprime,yprime)
        if posFunction[2][0]== np.Infinity:
            Ax+=0
        elif (cornerPos[1][1]<=cornerPos[2][1]):
            Ax += self.numeric_double_Integral(integrand[0],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],self.step,xprime,yprime)
        else:
            Ax += self.numeric_double_Integral(integrand[0],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[2][0]+posFunction[2][1],self.step,xprime,yprime)
        if posFunction[2][0]==np.Infinity:
            Ax+=0
        elif (posFunction[2][0]>=0):
            Ax += self.numeric_double_Integral(integrand[0],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],self.step,xprime,yprime)
        else:
            Ax += self.numeric_double_Integral(integrand[0],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[3][0]+posFunction[3][1],xprime*posFunction[2][0]+posFunction[2][1],self.step,xprime,yprime)
        if posFunction[0][0]== np.Infinity:
            Ay=0
        elif(posFunction[0][0]<=0):
            Ay = self.numeric_double_Integral(integrand[1],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],self.step,xprime,yprime)
        else:
            Ay = self.numeric_double_Integral(integrand[1],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[0][0]+posFunction[0][1],self.step,xprime,yprime)
        if posFunction[2][0]==np.Infinity:
            Ay+=0
        elif (cornerPos[1][1]<=cornerPos[2][1]):
            Ay += self.numeric_double_Integral(integrand[1],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],self.step,xprime,yprime)
        else:
            Ay += self.numeric_double_Integral(integrand[1],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[2][0]+posFunction[2][1],self.step,xprime,yprime)
        if posFunction[2][0]==np.Infinity:
            Ay+=0
        elif (posFunction[2][0]>=0):
            Ay += self.numeric_double_Integral(integrand[1],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],self.step,xprime,yprime)
        else:
            Ay += self.numeric_double_Integral(integrand[1],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[3][0]+posFunction[3][1],xprime*posFunction[2][0]+posFunction[2][1],self.step,xprime,yprime)
        if posFunction[0][0]== np.Infinity:
            Az=0
        elif(posFunction[0][0]<=0):
            Az = self.numeric_double_Integral(integrand[2],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],self.step,xprime,yprime)
        else:
            Az = self.numeric_double_Integral(integrand[2],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[0][0]+posFunction[0][1],self.step,xprime,yprime)
        if posFunction[2][0] == np.Infinity:
            Az+=0
        elif (cornerPos[1][1]<=cornerPos[2][1]):
            Az += self.numeric_double_Integral(integrand[2],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],self.step,xprime,yprime)
        else:
            Az += self.numeric_double_Integral(integrand[2],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[2][0]+posFunction[2][1],self.step,xprime,yprime)
        if posFunction[2][0]==np.Infinity:
            Az+=0
        elif (posFunction[2][0]>=0):
            Az += self.numeric_double_Integral(integrand[2],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],self.step,xprime,yprime)
        else:
            Az += self.numeric_double_Integral(integrand[2],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[3][0]+posFunction[3][1],xprime*posFunction[2][0]+posFunction[2][1],self.step,xprime,yprime)
        #The past series of conditional statements calculated magnetic potential, It does so by splitting the rectangular region into
        # 3 sections.

        Afunc = np.array([Ax,Ay,Az]) #creates vector of magnetic potetnial
        B = self.curl(Afunc,x,y,z)
        Bfun = lambdify([x,y,z],B)
        i=0 #loops evaluate anlytica B field at each point
        while i <int(2*self.xlim/self.resolution):
            j=0
            while j<int(2*self.ylim/self.resolution):
                k=0
                while k<int(2*self.zlim/self.resolution):

                    AijkX[i][j][k]=Bfun(i*self.resolution-self.xlim,j*self.resolution-self.ylim,k*self.resolution-self.zlim)[0]
                    AijkY[i][j][k]=Bfun(i*self.resolution-self.xlim,j*self.resolution-self.ylim,k*self.resolution-self.zlim)[1]
                    AijkZ[i][j][k]=Bfun(i*self.resolution-self.xlim,j*self.resolution-self.ylim,k*self.resolution-self.zlim)[2]

                    print("stuff",i,j,k,AijkX[i][j][k],AijkY[i][j][k],AijkZ[i][j][k])
                    k+=1
                j+=1
            i+=1

    # The 3d arrays holding the field values for each point

        threadQ.put(np.array([AijkX,AijkY,AijkZ])) # for handling multiple parts of the magnet
        return np.array([AijkX,AijkY,AijkZ])


    def Magnetic_field(self):

        AijkX=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        AijkY=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        AijkZ=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))

        threadlist=[]
        '''
        for i in range(int(2*self.diameter/(numofthreads*self.resolution))):
            j=i*numofthreads
            while j<numofthreads*(i+1):
                t=threading.Thread(target=self.Magnetic_Field_Slice,args=(np.sqrt( (self.diameter/2)**2-(self.resolution*j-self.diameter)**2 )  ,self.resolution*j-self.diameter))
                threadlist.append(t)
                t.start()
                j+=1
            while threadlist!=[]:
                threadlist.pop(0).join()
            while not threadQ.empty():
                A=threadQ.get()
                AijkX=np.add(AijkX,A[0])
                AijkY=np.add(AijkY,A[1])
                AijkZ=np.add(AijkZ,A[2])

        '''
        A=self.Magnetic_Field_Slice(self.diameter  ,0)
        AijkX=np.add(AijkX,A[0])
        AijkY=np.add(AijkY,A[1])
        AijkZ=np.add(AijkZ,A[2])
        return np.array([AijkX,AijkY,AijkZ])


    def sortbyx(self,positions):
        i=0
        while i<len(positions):
            j=len(positions)-i-1
            max = positions[j][0]
            maxi=j
            while j>=0:
                if positions[j][0]>max:
                    maxi=j
                    max= positions[j][0]
                j-=1
            placeholder = positions[len(positions)-i-1]
            positions[len(positions)-i-1]=positions[maxi]
            positions[maxi]= placeholder
            i+=1
        return positions
        # function sorts a vector list based on the x componet of a vector
        #used for the currentPos function inorder to corecly map postions in all regions
    def currentPos(self,thickness,diameter,armLength,phi): #finds the position of each of the four cornes of the magnet
        psi =  np.arctan(diameter/(2*(armLength)))    #angle made between center of magnet and top front end
        theta = np.arctan(diameter/(2*(armLength+thickness))) # angle made between center of the magnet and backend

        alpha = (armLength**2 + (diameter/2)**2)**.5 # distance from front end
        beta = ((armLength+thickness)**2 + (diameter/2)**2)**.5 # distance from back end
        r0= np.array( [ alpha*np.cos(phi+psi),alpha*np.sin(phi+psi),0] )
        r1= np.array( [ alpha*np.cos(phi-psi),alpha*np.sin(phi-psi),0] )
        r2= np.array( [ beta*np.cos(phi+theta),beta*np.sin(phi+theta),0] )
        r3= np.array( [ beta*np.cos(phi-theta),beta*np.sin(phi-theta),0] )
        return [r0,r1,r2,r3]
    #return sortbyx([r0,r1,r2,r3])
    def postionFunctions(self,thickness,diameter,armLength,phi):
        L = self.sortbyx(self.currentPos(thickness,diameter,armLength,phi))

        tolerance =500
        if (L[0][0]-L[2][0])!=0:
            slopeofLine0To2 = (L[0][1]-L[2][1])/(L[0][0]-L[2][0])
            if slopeofLine0To2>tolerance:
                slopeofLine0To2=np.Infinity

        else:
            slopeofLine0To2 = np.Infinity
        if (L[0][0]-L[1][0])!=0:
            slopeofLine0To1 = (L[0][1]-L[1][1])/(L[0][0]-L[1][0])
            if slopeofLine0To1>tolerance:
                slopeofLine0To1=np.Infinity
        else:
            slopeofLine0To1 = np.Infinity
        if (L[1][0]-L[3][0])!=0:
            slopeofLine1To3 = (L[1][1]-L[3][1])/(L[1][0]-L[3][0])
            if slopeofLine1To3>tolerance:
                slopeofLine1To3=np.Infinity
        else:
            slopeofLine1To3 = np.Infinity
        if (L[2][0]-L[3][0])!=0:
            slopeofLine2To3 =  (L[2][1]-L[3][1])/(L[2][0]-L[3][0])
            if slopeofLine2To3>tolerance:
                slopeofLine2To3=np.Infinity
        else:
            slopeofLine2To3 = np.Infinity

        b02 = L[0][1]-slopeofLine0To2*L[0][0]
        b01 = L[0][1]-slopeofLine0To1*L[0][0]
        b13 = L[1][1]-slopeofLine1To3*L[1][0]
        b23 = L[2][1]-slopeofLine2To3*L[2][0]
    #print((round(L[0][0],2),round(L[0][1],2) ),(round(L[1][0],2),round(L[1][1],2) ),(round(L[2][0],2),round(L[2][1],2) ),(round(L[3][0],2),round(L[3][1],2) ) )
    #print("current pos",(round(slopeofLine0To1,2),round(b01,2)),(round(slopeofLine0To2,2),round(b02,2)),(round(slopeofLine1To3,2),round(b13,2)),(round(slopeofLine2To3,2),round(b23,2)))
        return [(slopeofLine0To1,b01),(slopeofLine0To2,b02),(slopeofLine1To3,b13),(slopeofLine2To3,b23)]
    def numeric_double_Integral(self,integrand,x0,xf,yb,yt,xstep,ystep,symx,symy):

    #simple rieman double integral
        xi=x0
        I=0.0
        f= lambdify(symx,yb,"numpy")
        g= lambdify(symx,yt,"numpy")
        h = lambdify([symx,symy],integrand)
        while( xi<xf):
            yi=f(xi)
            while( yi< g(xi)):
                area = h(xi+xstep*.5,yi+ystep*.5)
                I+=area*ystep*xstep

                yi+=ystep
            xi+=xstep
        return I
    def numeric_curl(Vx,Vy,Vz,step): #legacy function, used to calculate curl
        #numeric curl function, using linear approximation
        Wx=np.zeros(Vx.shape)
        Wy=np.zeros(Vy.shape)
        Wz=np.zeros(Vz.shape)
        print(Vx.shape)
        print(Vy.shape)
        print(Vz.shape)
        for i in range(Vx.shape[0]):
            for j in range(Vy.shape[0]):
                for k in range(Vz.shape[0]):
                    if j+1<Vz.shape[0] and k+1<Vz.shape[0]:
                        Wx[i][j][k]=(Vz[i][j+1][k]-Vz[i][j][k])/step+(Vy[i][j][k+1]-Vy[i][j][k])/step
                    elif k+1<Vz.shape[0]:
                        Wx[i][j][k]= (Vy[i][j][k+1]-Vy[i][j][k])/step
                    elif j+1<Vz.shape[0]:
                        Wx[i][j][k]= (Vz[i][j+1][k]-Vz[i][j][k])/step
                    else:
                         Wy[i][j][k]=Wy[i][j-1][k-1]
                    if i+1<Vz.shape[0] and k+1<Vz.shape[0]:
                        Wy[i][j][k]=-1*((Vz[i+1][j][k]-Vz[i][j][k])/step+(Vx[i][j][k+1]-Vx[i][j][k])/step)
                    elif k+1<Vz.shape[0]:
                        Wy[i][j][k] = ((Vx[i][j][k+1]-Vx[i][j][k])/step)
                    elif i+1<Vz.shape[0]:
                        Wy[i][j][k]= -1*((Vz[i+1][j][k]-Vz[i][j][k])/step)
                    else:
                         Wy[i][j][k]=Wy[i-1][j][k-1]
                    if j+1<Vz.shape[0] and i+1<Vz.shape[0]:
                        Wz[i][j][k]= (Vy[i+1][j][k]-Vy[i][j][k])/step+(Vx[i][j+1][k]-Vx[i][j][k])/step
                    elif j+1<Vz.shape[0]:
                        Wz[i][j][k]= (Vx[i][j+1][k]-Vx[i][j][k])/step
                    elif i+1<Vz.shape[0]:
                        Wz[i][j][k]= (Vy[i+1][j][k]-Vy[i][j][k])/step
                    else:
                         Wz[i][j][k]=Wz[i-1][j-1][k]
        return np.array([Wx,Wy,Wz])



    def rotating_field(self):
        #total time is the period of time, this program takes time slices of rotating magnets
        numberofticks= int(self.totaltime/self.timestep)
        Bfield=[] # bfield
        for t in range(numberofticks):
            print("started",t,"!")
            start = time.time()
            self.phi=self.angularvelocity*t+self.startphi
            B=self.Magnetic_field()
            #Bfield.append(Magnetic_Field(A,step))
            file = open("Bfielddata"+str(t)+".txt","w+")

            for i in range(np.shape(B[0])[0]):
                for j in range(np.shape(B[0])[0]):
                    for k in range(np.shape(B[0])[0]):
                        p= str(i*self.resolution-self.xlim)
                        q= str(j*self.resolution-self.ylim)
                        r= str(k*self.resolution-self.zlim)
                        a=str(B[0][i][j][k])
                        b=str(B[1][i][j][k])
                        c=str(B[2][i][j][k])
                        print("("+p+","+q+","+r+")  ("+a+","+b+","+c+")")
                        file.write("("+p+","+q+","+r+")  ("+a+","+b+","+c+")\r\n")
            delt=time.time()-start
            print("Finished ",t,"! It took",delt," seconds")

        return Bfield
    def numeric_double_Integral(self,integrand,x0,xf,yb,yt,step,symx,symy):
        #print("Entering x's", x0,xf)
        xstep=step
        ystep=step
        xi=x0
        I=0
        #print("funcs",yb,yt)
        f= lambdify(symx,yb,"numpy")
        g= lambdify(symx,yt,"numpy")
        h = lambdify([symx,symy],integrand)
        while( round(xi,6)<round(xf,6)):
            yi=f(xi)
            #print("x,yb and yt",x,y,g(x),f(x))
            while( round(yi,6)< round(g(xi),6)):

                area = h(xi+xstep*.5,yi+ystep*.5)

                I+=area*ystep*xstep
                #print("I,f,g,", I,y,g(x))
                yi+=ystep
            xi+=xstep
        #print("Returned I,", I)
        return I

res =.5
xl =12
yl =12
zl =12
stp=.25
tt=48.0
ts=1.0
sphi =0*np.pi/5
th=2
dia=4
aL=3
av=np.pi/24
m=1.0
M = MARS(th,dia,aL,sphi, av,m,stp,res ,xl, yl,zl, tt,ts)
M.rotating_field()
