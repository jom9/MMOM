from sympy import *
import numpy as np
import time
from queue import Queue
import threading
import sys
from scipy import constants
numofthreads=40
threadQ=Queue()

def magnitude(vector):
    m=0
    for i in range(vector.size):
        m+=vector[i]**2
    return m**.5
def eval_field(A,symx,symy,symz):
    B=[]
    B+=[diff(A[2],symy)-diff(A[1],symz)]
    B+=[-1*(diff(A[2],symx)-diff(A[0],symz))]
    B+=[diff(A[1],symx)-diff(A[0],symy)]

    return B
def Magnetic_Potenial(thickness,diameter,armLength,phi, magnetization,step,zprime,xlim, ylim,zlim,resolution):
    # thickness is the thickness of the magnet
    # diameter of the magnet
    #armLength is the distance the magnet is away from the center
    #phi is the current angular position
    #magnetization is the magnatic dipole moments per unit area, this is assumed to be coaxial
    #returns the analyitcal function of the magentic potential
    M = magnetization*np.array([np.cos(phi),np.sin(phi),0]) # magnetization vector

    x,y,z,xprime,yprime= symbols('x y z xprime yprime')
    delR = np.array([x-xprime,y-yprime,z-zprime]) # this is the distance vector, with the primed values denoting source points and the unprimed values denoting field values

    integrand =constants.mu_o/(4*np.pi)np.cross(M,delR)/(magnitude(delR)**3)
    cornerPos= sortbyx(currentPos(thickness,diameter,armLength,phi)) #gets the four corners needed to take each rectangluar slice
    posFunction = postionFunctions(thickness,diameter,armLength,phi) #creates the lines that trace out borders to slice

    AijkX=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    AijkY=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    AijkZ=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))

    if posFunction[0][0]== np.Infinity:
        Ax =0
    elif(posFunction[0][0]<=0):
        Ax = numeric_double_Integral(integrand[0],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    else:
        Ax = numeric_double_Integral(integrand[0],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[0][0]+posFunction[0][1],step,xprime,yprime)
    if posFunction[2][0]== np.Infinity:
        Ax+=0
    elif (cornerPos[1][1]<=cornerPos[2][1]):
        Ax += numeric_double_Integral(integrand[0],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    else:
        Ax += numeric_double_Integral(integrand[0],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
    if posFunction[2][0]==np.Infinity:
        Ax+=0
    elif (posFunction[2][0]>=0):
        Ax += numeric_double_Integral(integrand[0],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)
    else:
        Ax += numeric_double_Integral(integrand[0],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[3][0]+posFunction[3][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
    if posFunction[0][0]== np.Infinity:
        Ay=0
    elif(posFunction[0][0]<=0):
        Ay = numeric_double_Integral(integrand[1],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    else:
        Ay = numeric_double_Integral(integrand[1],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[0][0]+posFunction[0][1],step,xprime,yprime)
    if posFunction[2][0]==np.Infinity:
        Ay+=0
    elif (cornerPos[1][1]<=cornerPos[2][1]):
        Ay += numeric_double_Integral(integrand[1],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    else:
        Ay += numeric_double_Integral(integrand[1],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
    if posFunction[2][0]==np.Infinity:
        Ay+=0
    elif (posFunction[2][0]>=0):
        Ay += numeric_double_Integral(integrand[1],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)
    else:
        Ay += numeric_double_Integral(integrand[1],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[3][0]+posFunction[3][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
    if posFunction[0][0]== np.Infinity:
        Az=0
    elif(posFunction[0][0]<=0):
        Az = numeric_double_Integral(integrand[2],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    else:
        Az = numeric_double_Integral(integrand[2],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[0][0]+posFunction[0][1],step,xprime,yprime)
    if posFunction[2][0] == np.Infinity:
        Az+=0
    elif (cornerPos[1][1]<=cornerPos[2][1]):
        Az += numeric_double_Integral(integrand[2],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    else:
        Az += numeric_double_Integral(integrand[2],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
    if posFunction[2][0]==np.Infinity:
        Az+=0
    elif (posFunction[2][0]>=0):
        Az += numeric_double_Integral(integrand[2],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)
    else:
        Az += numeric_double_Integral(integrand[2],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[3][0]+posFunction[3][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
    # evualtes the integral to get the potential. each component takes 3 integrals since the region must be broken up into 3 subregions

    Afunc = np.array([Ax,Ay,Az])
    B = eval_field(Afunc,x,y,z)
    Bfun = lambdify([x,y,z],B)
    i=0
    while i <int(2*xlim/resolution):
        j=0
        while j<int(2*ylim/resolution):
            k=0
            while k<int(2*zlim/resolution):

                AijkX[i][j][k]=Bfun(i*resolution-xlim,j*resolution-ylim,k*resolution-zlim)[0]
                AijkY[i][j][k]=Bfun(i*resolution-xlim,j*resolution-ylim,k*resolution-zlim)[1]
                AijkZ[i][j][k]=Bfun(i*resolution-xlim,j*resolution-ylim,k*resolution-zlim)[2]
                print("stuff",i,j,k,AijkX[i][j][k],AijkY[i][j][k],AijkZ[i][j][k])
                k+=1
            j+=1
        i+=1

    # The 3d arrays holding the field values for each point

    threadQ.put(np.array([AijkX,AijkY,AijkZ]))
    return np.array([AijkX,AijkY,AijkZ])


def Magnetic_Potenial_field(thickness,diameter,armLength,phi, magnetization,step,resolution ,xlim, ylim,zlim):
    k= zlim*-1
    AijkX=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    AijkY=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    AijkZ=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))

    threadlist=[]
    '''
    for i in range(int(2*zlim/(numofthreads*resolution))):
        j=i*numofthreads
        while j<numofthreads*(i+1):
            t=threading.Thread(target=Magnetic_Potenial,args=(thickness,diameter,armLength,phi, magnetization,step,resolution*j-zlim,xlim, ylim,zlim,resolution))
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
    A=Magnetic_Potenial(thickness,diameter,armLength,phi, magnetization,step,0,xlim, ylim,zlim,resolution)
    AijkX=np.add(AijkX,A[0])
    AijkY=np.add(AijkY,A[1])
    AijkZ=np.add(AijkZ,A[2])
    return np.array([AijkX,AijkY,AijkZ])

    '''
    #this generates the entire potential over the entire magnets region by dividing it into slices go from the negitve z limit to the positive z limit
    k= zlim*-1
    throwaway= symbols('throwaway')
    M=np.array([throwaway,throwaway,throwaway],dtype='object')
    M-=throwaway


    while k<zlim:
        Aith=Magnetic_Potenial(thickness,diameter,armLength,phi, magnetization,step,k)
        M[0]+=Aith[0]
        M[1]+=Aith[1]
        M[2]+=Aith[2]
        k   +=step
    x,y,z =symbols('x y z')
    A= lambdify([x,y,z],M) #analyitcal potential function geneated by combining slices

    AijkX=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    AijkY=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    AijkZ=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    # The 3d arrays holding the field values for each point
    i=0
    while i <int(2*xlim/resolution):
        j=0
        while j<int(2*ylim/resolution):
            k=0
            while k<int(2*zlim/resolution):

                AijkX[i][j][k]=A(i*resolution-xlim,j*resolution-ylim,k*resolution-zlim)[0]
                AijkY[i][j][k]=A(i*resolution-xlim,j*resolution-ylim,k*resolution-zlim)[1]
                AijkZ[i][j][k]=A(i*resolution-xlim,j*resolution-ylim,k*resolution-zlim)[2]
                print("stuff",i,j,k,AijkX[i][j][k],AijkY[i][j][k],AijkZ[i][j][k])
                k+=1
            j+=1
        i+=1
    return np.array([AijkX,AijkY,AijkZ])
    '''
def sortbyx(positions):
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


def currentPos(thickness,diameter,armLength,phi): #finds the position of each of the four cornes of the magnet
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
def postionFunctions(thickness,diameter,armLength,phi):
    L = sortbyx(currentPos(thickness,diameter,armLength,phi))

    tolerance = 500
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
def numeric_double_Integral(integrand,x0,xf,yb,yt,xstep,ystep,symx,symy):

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
def numeric_curl(Vx,Vy,Vz,step):
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

def Magnetic_Field(potential,step):
    B=numeric_curl(potential[0],potential[1],potential[2],step)
    return B
def double_integral_checker(thickness,diameter,armLength,startphi, angularvelocity,magnetization,step,resolution ,xlim, ylim,zlim, totaltime,timestep):
        tolerance= 500
        numberofticks= int(totaltime/timestep)
        Bfield=[] # bfield
        zprime =0
        for t in range(numberofticks):
            print("started",t,"!")
            start = time.time()
            phi=angularvelocity*t+startphi

            M = magnetization*np.array([np.cos(phi),np.sin(phi),0]) # magnetization vector

            x,y,z,xprime,yprime= symbols('x y z xprime yprime')
            delR = np.array([x-xprime,y-yprime,z-zprime]) # this is the distance vector, with the primed values denoting source points and the unprimed values denoting field values

            integrand =np.cross(M,delR)/(magnitude(delR)**3)

            cornerPos= sortbyx(currentPos(thickness,diameter,armLength,phi)) #gets the four corners needed to take each rectangluar slice
            posFunction = postionFunctions(thickness,diameter,armLength,phi) #creates the lines that trace out borders to slice


            if posFunction[0][0]== np.Infinity:
                Ax =0
            elif(posFunction[0][0]<=0):
                Ax = numeric_double_Integral(integrand[0],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            else:
                Ax = numeric_double_Integral(integrand[0],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[0][0]+posFunction[0][1],step,xprime,yprime)
            if posFunction[2][0]== np.Infinity:
                Ax+=0
            elif (cornerPos[1][1]<=cornerPos[2][1]):
                Ax += numeric_double_Integral(integrand[0],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            else:
                Ax += numeric_double_Integral(integrand[0],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
            if posFunction[2][0]==np.Infinity:
                Ax+=0
            elif (posFunction[2][0]>=0):
                Ax += numeric_double_Integral(integrand[0],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)
            else:
                Ax += numeric_double_Integral(integrand[0],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[3][0]+posFunction[3][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
            if posFunction[0][0]== np.Infinity:
                Ay=0
            elif(posFunction[0][0]<=0):
                Ay = numeric_double_Integral(integrand[1],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            else:
                Ay = numeric_double_Integral(integrand[1],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[0][0]+posFunction[0][1],step,xprime,yprime)
            if posFunction[2][0]==np.Infinity:
                Ay+=0
            elif (cornerPos[1][1]<=cornerPos[2][1]):
                Ay += numeric_double_Integral(integrand[1],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            else:
                Ay += numeric_double_Integral(integrand[1],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
            if posFunction[2][0]==np.Infinity:
                Ay+=0
            elif (posFunction[2][0]>=0):
                Ay += numeric_double_Integral(integrand[1],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)
            else:
                Ay += numeric_double_Integral(integrand[1],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[3][0]+posFunction[3][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
            if posFunction[0][0]== np.Infinity:
                Az=0
            elif(posFunction[0][0]<=0):
                Az = numeric_double_Integral(integrand[2],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            else:
                Az = numeric_double_Integral(integrand[2],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[0][0]+posFunction[0][1],step,xprime,yprime)
            if posFunction[2][0] == np.Infinity:
                Az+=0
            elif (cornerPos[1][1]<=cornerPos[2][1]):
                Az += numeric_double_Integral(integrand[2],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            else:
                Az += numeric_double_Integral(integrand[2],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[1][0]+posFunction[1][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
            if posFunction[2][0]==np.Infinity:
                Az+=0
            elif (posFunction[2][0]>=0):
                Az += numeric_double_Integral(integrand[2],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)
            else:
                Az += numeric_double_Integral(integrand[2],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[3][0]+posFunction[3][1],xprime*posFunction[2][0]+posFunction[2][1],step,xprime,yprime)
            # evualtes the integral to get the potential. each component takes 3 integrals since the region must be broken up into 3 subregions
            A=np.array([Ax,Ay,Az])
            print(phi,A)
            delt=time.time()-start
            print("Finished ",t,"! It took",delt," seconds")

        return Bfield
def rotating_field(thickness,diameter,armLength,startphi, angularvelocity,magnetization,step,resolution ,xlim, ylim,zlim, totaltime,timestep):
        #total time is the period of time, this program takes time slices of rotating magnets
        numberofticks= int(totaltime/timestep)
        Bfield=[] # bfield
        for t in range(numberofticks):
            print("started",t,"!")
            start = time.time()
            phi=angularvelocity*t+startphi
            B=Magnetic_Potenial_field(thickness,diameter,armLength,phi, magnetization,step,resolution ,xlim, ylim,zlim)
            #Bfield.append(Magnetic_Field(A,step))
            file = open("Bfielddata"+str(t)+".txt","w+")

            for i in range(np.shape(B[0])[0]):
                for j in range(np.shape(B[0])[0]):
                    for k in range(np.shape(B[0])[0]):
                        p= str(i*res-xlim)
                        q= str(j*res-ylim)
                        r= str(k*res-zlim)
                        a=str(B[0][i][j][k])
                        b=str(B[1][i][j][k])
                        c=str(B[2][i][j][k])
                        print("("+p+","+q+","+r+")  ("+a+","+b+","+c+")")
                        file.write("("+p+","+q+","+r+")  ("+a+","+b+","+c+")\r\n")
            delt=time.time()-start
            print("Finished ",t,"! It took",delt," seconds")

        return Bfield

def numeric_double_Integral(integrand,x0,xf,yb,yt,step,symx,symy):
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
xlim =12
ylim =12
zlim =12
step=.25
totaltime=48.0
timestep=1.0
startphi =0*np.pi/5
thickness=2
diameter=4
armLength=3
angularvelocity=np.pi/24
magnetization=1.0
#rotating_field(thickness,diameter,armLength,startphi, angularvelocity,magnetization,step,resolution ,xlim, ylim,zlim, totaltime,timestep):
#double_integral_checker(thickness,diameter,armLength,startphi, angularvelocity,magnetization,step,res ,xlim, ylim,zlim, totaltime,timestep)
B =rotating_field(thickness,diameter,armLength,startphi, angularvelocity,magnetization,step,res ,xlim, ylim,zlim, totaltime,timestep)
'''
for l in range(B.size):
    file = open("Bfielddata"+str(l)+".txt","w+")
    A=B[l]
    for i in range(np.shape(A[0])[0]):
        for j in range(np.shape(A[0])[0]):
            for k in range(np.shape(A[0])[0]):
                p= str(i*res-xlim)
                q= str(j*res-ylim)
                r= str(k*res-zlim)
                a=str(A[0][i][j][k])
                b=str(A[1][i][j][k])
                c=str(A[2][i][j][k])
                print("("+p+","+q+","+r+")  ("+a+","+b+","+c+")")
                file.write("("+p+","+q+","+r+")  ("+a+","+b+","+c+")\r\\n")
'''
