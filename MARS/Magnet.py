import numpy as np
from sympy import *
from scipy import constants
import os
class Magnet():
    def __init__(self,thickness,diameter,armLength,magnetization,xlim,ylim,zlim,resolution,step):
        self.thickness=thickness
        self.diameter=diameter
        self.armLength=armLength
        self.xlim=xlim
        self.ylim=ylim
        self.zlim=zlim
        self.resolution=resolution
        self.step=step
        self.magnetization=magnetization
    def Magnetic_Field_Slice(self,slice_diameter,zprime,phi):
    # thickness is the thickness of the magnet
    # diameter of the magnet
    #armLength is the distance the magnet is away from the center
    #phi is the current angular position
    #magnetization is the magnatic dipole moments per unit area, this is assumed to be coaxial
    #returns the 3d array of the magentic field

        M = self.magnetization*np.array([np.cos(phi),np.sin(phi),0]) # magnetization vector

        x,y,z,xprime,yprime= symbols('x y z xprime yprime')
        delR = np.array([x-xprime,y-yprime,z-zprime]) # this is the distance vector, with the primed values denoting source points and the unprimed values denoting field values

        integrand =np.cross(M,delR)/(self.magnitude(delR)**3)

        cornerPos= self.sortbyx(self.currentPos(self.thickness,slice_diameter,self.armLength,phi)) #gets the four corners needed to take each rectangluar slice

        posFunction = self.postionFunctions(self.thickness,slice_diameter,self.armLength,phi) #creates the lines that trace out borders to slice

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

        #threadQ.put(np.array([AijkX,AijkY,AijkZ])) # for handling multiple parts of the magnet
        return np.array([AijkX,AijkY,AijkZ])



    def genBfield(self,phi):

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
        A=self.Magnetic_Field_Slice(self.diameter  ,0,phi)
        AijkX=np.add(AijkX,A[0])
        AijkY=np.add(AijkY,A[1])
        AijkZ=np.add(AijkZ,A[2])

        return np.array([AijkX,AijkY,AijkZ])
    def genBfield(self,phi,i):
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
        A=self.Magnetic_Field_Slice(self.diameter  ,0,phi)
        AijkX=np.add(AijkX,A[0])
        AijkY=np.add(AijkY,A[1])
        AijkZ=np.add(AijkZ,A[2])
        current_path = os.getcwd()
        if "magnet_Bfielddata" not in os.listdir(current_path):
            os.mkdir("magnet_Bfielddata")
        os.chdir("magnet_Bfielddata")

        file = open("magnet_Bfielddata"+str(i)+".txt","w+")

        for i in range(int(2*self.xlim/self.resolution)):
            for j in range(int(2*self.ylim/self.resolution)):
                for k in range(int(2*self.zlim/self.resolution)):
                    file.write( "("+str(i*self.resolution-self.xlim)+","+str(j*self.resolution-self.ylim)+","+str(k*self.resolution-self.zlim)+")  "+"("+str(AijkX[i][j][k])+","+str(AijkY[i][j][k])+","+str(AijkZ[i][j][k]) +")\r\n")
        os.chdir(current_path)
        return np.array([AijkX,AijkY,AijkZ])
    def getBfield(self,i):

        try:
            current_path=os.getcwd()
            os.chdir("magnet_Bfielddata")
            file = open("magnet_Bfielddata"+str(i)+".txt",'r')
            Bx_field = np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
            By_field = np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
            Bz_field = np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
            for line in file.readlines():
                parsedline = line.split(" ")

                coors = parsedline[0].split(",")
                field = parsedline[-1].split(",")
                if len(coors)==3 and len(field)==3:

                    print(i,coors,field)
                    x = float(coors[0][1:])
                    y = float(coors[1])
                    z = float(coors[2][:-1])

                    u = float(field[0][1:])
                    v = float(field[1])
                    w = float(field[2][:-2])
                    Bx_field[int(x//self.resolution-self.xlim)][int(y//self.resolution-self.ylim)][int(z//self.resolution-self.zlim)]=u
                    By_field[int(x//self.resolution-self.xlim)][int(y//self.resolution-self.ylim)][int(z//self.resolution-self.zlim)]=v
                    Bz_field[int(x//self.resolution-self.xlim)][int(y//self.resolution-self.ylim)][int(z//self.resolution-self.zlim)]=w
                else:
                    continue
            os.chdir(current_path)
        except NotADirectoryError:
            print("Nothing to get")

        return np.array([Bx_field,By_field,Bz_field])
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
    def isInMagnet(self,x,y,phi):
        epsilon = 6e-17
        r = self.currentPos(self.thickness,self.diameter,self.armLength,phi)
        dist_point_0 = ((((r[0][0] - x) ** 2) + (((r[0][1] - y)) ** 2))) ** (0.5)
        dist_point_1 = ((((r[1][0] - x) ** 2) + (((r[1][1] - y)) ** 2))) ** (0.5)
        dist_point_2 = ((((r[2][0] - x) ** 2) + (((r[2][1] - y)) ** 2))) ** (0.5)
        dist_point_3 = ((((r[3][0] - x) ** 2) + (((r[3][1] - y)) ** 2))) ** (0.5)
        semiperimeterp02 = ((self.thickness)+(dist_point_0)+(dist_point_2)) / 2
        semiperimeterp01 = ((self.diameter)+(dist_point_0)+(dist_point_1)) / 2
        semiperimeterp23 = ((self.diameter)+(dist_point_2)+(dist_point_3)) / 2
        semiperimeterp13 = ((self.thickness)+(dist_point_1)+(dist_point_3)) / 2
        Areap02 = ((semiperimeterp02)*((semiperimeterp02) - (dist_point_0))*((semiperimeterp02) - (dist_point_2))*((semiperimeterp02) - (self.thickness)))**(0.5)
        Areap01 = ((semiperimeterp01)*((semiperimeterp01) - (dist_point_0))*((semiperimeterp01) - (dist_point_1))*((semiperimeterp01) - (self.diameter)))**(0.5)
        Areap23 = ((semiperimeterp23)*((semiperimeterp23) - (dist_point_2))*((semiperimeterp23) - (dist_point_3))*((semiperimeterp23) - (self.diameter)))**(0.5)
        Areap13 = ((semiperimeterp13)*((semiperimeterp13) - (dist_point_1))*((semiperimeterp13) - (dist_point_3))*((semiperimeterp13) - (self.thickness)))**(0.5)
        SumArea = Areap01 + Areap02 + Areap23 + Areap13 + epsilon
        if SumArea <= Areaofmagnet:
            return True
        else:
            return False
    def genForce(self,B_field,phi,i):
        Bx_field=B_field[0]
        By_field=B_field[1]
        Bz_field=B_field[2]
        current_path = os.getcwd()
        if "forcefielddata" not in os.listdir(current_path):
            os.mkdir("forcefielddata")
        os.chdir("forcefielddata")
        file = open("forcefielddata"+str(i)+".txt","w+")
        forcex=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        forcey=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        forcez=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        M = self.magnetization*np.array([np.cos(phi),np.sin(phi),0])
        for i in range(int(2*self.xlim/self.resolution)):
            for j in range(int(2*self.ylim/self.resolution)):
                for k in range(int(2*self.zlim/self.resolution)):
                    force= np.grad( np.dot(M,np.array([Bx_field[i][j][k],By_field[i][j][k],Bz_field[i][j][k]  ])))
                    forcex[i][j][k]=force[0]
                    forcey[i][j][k]=force[1]
                    forcez[i][j][k]=force[2]
                    file.write( "("+str(i*self.resolution-self.xlim)+","+str(j*self.resolution-self.ylim)+","+str(k*self.resolution-self.zlim)+")  "+"("+str(force[0])+","+str(force[1])+","+str(force[2]) +")\r\n")
        os.chdir(current_path)
    def getForce(self,i):
        forcex=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        forcey=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        forcez=np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
        current_path = os.getcwd()
        if "forcefielddata" not in os.listdir(current_path):
            os.mkdir("forcefielddata")
        os.chdir("forcefielddata")
        file = open("forcefielddata"+str(i)+".txt","r")
        for line in file.readlines():
            parsedline = line.split(" ")

            coors = parsedline[0].split(",")
            field = parsedline[-1].split(",")
            if len(coors)==3 and len(field)==3:

                print(i,coors,field)
                x = float(coors[0][1:])
                y = float(coors[1])
                z = float(coors[2][:-1])

                u = float(field[0][1:])
                v = float(field[1])
                w = float(field[2][:-2])
                forcex[int(x//self.resolution-self.xlim)][int(y//self.resolution-self.ylim)][int(z//self.resolution-self.zlim)]=u
                forcey[int(x//self.resolution-self.xlim)][int(y//self.resolution-self.ylim)][int(z//self.resolution-self.zlim)]=v
                forcez[int(x//self.resolution-self.xlim)][int(y//self.resolution-self.ylim)][int(z//self.resolution-self.zlim)]=w
            else:
                continue
        os.chdir(current_path)
        return np.array([forcex,forcey,forcez])
