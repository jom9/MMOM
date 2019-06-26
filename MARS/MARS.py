from sympy import *
import numpy as np
def magnitude(vector):
    m=0
    for i in range(vector.size):
        m+=vector[i]**2
    return m**.5
def Magnetic_Potenial(thickness,diameter,armLength,phi, magnetization,step,zprime):
    # thickness is the thickness of the magnet
    # diameter of the magnet
    #armLength is the distance the magnet is away from the center
    #phi is the current angular position
    #magnetization is the magnatic dipole moments per unit area, this is assumed to be coaxial
    #returns the analyitcal function of the magentic potential
    M = np.array([np.cos(phi),np.sin(phi),0]) # magnetization vector

    x,y,z,xprime,yprime= symbols('x y z xprime yprime')
    delR = np.array([x-xprime,y-yprime,z-zprime]) # this is the distance vector, with the primed values denoting source points and the unprimed values denoting field values

    integrand = np.cross(M,delR)/(magnitude(delR)**3)
    
    cornerPos= currentPos(thickness,diameter,armLength,phi)
    posFunction = postionFunctions(thickness,diameter,armLength,phi)

    Ax = numeric_double_Integral(integrand[0],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Ax += numeric_double_Integral(integrand[0],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Ax += numeric_double_Integral(integrand[0],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)

    Ay = numeric_double_Integral(integrand[1],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Ay += numeric_double_Integral(integrand[1],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Ay += numeric_double_Integral(integrand[1],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)

    Az = numeric_double_Integral(integrand[2],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Az += numeric_double_Integral(integrand[2],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Az += numeric_double_Integral(integrand[2],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)

    A=np.array([Ax,Ay,Az])
    return np.array([Ax,Ay,Az])

def Magnetic_Potenial_field(thickness,diameter,armLength,phi, magnetization,step,resolution ,xlim, ylim,zlim):
    M=np.array([0,0,0])
    k= zlim*-1
    while k<zlim:
        M+=Magnetic_Potenial(thickness,diameter,armLength,phi, magnetization,step,k)
        k+=step
    x,y,z =symbols('x y z')
    A= lambdify([x,y,z],M)

    AijX=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    AijY=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    AijZ=np.zeros((int(2*xlim/resolution),int(2*ylim/resolution),int(2*zlim/resolution)))
    i=0
    while i <int(2*xlim/resolution):
        j=0
        while j<int(2*ylim/resolution):
            k=0
            while k<int(2*zlim/resolution):
                AijX[i][j][k]=A(i*resolution-xlim,j*resolution-ylim)[0]
                AijY[i][j][k]=A(i*resolution-xlim,j*resolution-ylim)[1]
                AijZ[i][j][k]=A(i*resolution-xlim,j*resolution-ylim)[2]
                k+=1
            j+=1
        i+=1
    return np.array([AijX,AijY,AijZ])

def currentPos(thickness,diameter,armLength,phi): #finds the position of each of the four cornes of the magnet
    psi = np.absolute( np.arctan(diameter/2*(armLength))   ) #angle made between center of magnet and top front end
    theta = np.absolute( np.arctan(diameter/2*(armLength+thickness))) # angle made between center of the magnet and backend

    alpha = (armLength**2 + (diameter/2)**2)**1/2 # distance from front end
    beta = ((armLength+thickness)**2 + (diameter/2)**2)**1/2 # distance from back end
    r0= np.array( [ alpha*np.cos(phi+psi),alpha*np.sin(phi+psi),0] )
    r1= np.array( [ alpha*np.cos(phi-psi),alpha*np.sin(phi-psi),0] )
    r2= np.array( [ beta*np.cos(phi+theta),beta*np.sin(phi+theta),0] )
    r3= np.array( [ beta*np.cos(phi-theta),beta*np.sin(phi-theta),0] )
    return [r0,r1,r2,r3]
def postionFunctions(thickness,diameter,armLength,phi):
    L = currentPos(thickness,diameter,armLength,phi)
    if (L[0][0]-L[2][0])!=0:
        slopeofLine0To2 = (L[0][1]-L[2][1])/(L[0][0]-L[2][0])
    else:
        slopeofLine0To2 = np.Infinity
    if (L[0][0]-L[1][0])!=0:
        slopeofLine0To1 = (L[0][1]-L[1][1])/(L[0][0]-L[1][0])
    else:
        slopeofLine0To1 = np.Infinity
    if (L[1][0]-L[3][0])!=0:
        slopeofLine1To3 = (L[1][1]-L[3][1])/(L[1][0]-L[3][0])
    else:
        slopeofLine1To3 = np.Infinity
    if (L[2][0]-L[3][0])!=0:
        slopeofLine2To3 =  (L[2][1]-L[3][1])/(L[2][0]-L[3][0])
    else:
        slopeofLine2To3 = np.Infinity

    b02 = L[0][1]-slopeofLine0To2*L[0][0]
    b01 = L[0][1]-slopeofLine0To1*L[0][0]
    b13 = L[1][1]-slopeofLine1To3*L[1][0]
    b23 = L[2][1]-slopeofLine2To3*L[2][0]


    return [(slopeofLine0To1,b01),(slopeofLine0To2,b02),(slopeofLine1To3,b13),(slopeofLine2To3,b23)]
def numeric_double_Integral(integrand,x0,xf,yb,yt,xstep,ystep,symx,symy):
    x=x0
    I=0
    f= lambdify(symx,yb,"numpy")
    g= lambdify(symx,yt,"numpy")
    h = lambdify([symx,symy],integrand)
    while( x<xf):
        y=f(x)
        while( y< g(x)):
            area = h(x+xstep*.5,y+ystep*.5)

            I+=area*ystep*xstep

            y+=ystep
        x+=xstep
    return I

def numeric_double_Integral(integrand,x0,xf,yb,yt,step,symx,symy):
    xstep=step
    ystep=step
    x=x0
    I=0
    f= lambdify(symx,yb,"numpy")
    g= lambdify(symx,yt,"numpy")
    h = lambdify([symx,symy],integrand)
    while( x<xf):
        y=f(x)
        while( y< g(x)):
            area = h(x+xstep*.5,y+ystep*.5)

            I+=area*ystep*xstep

            y+=ystep
        x+=xstep
    return I

def numeric_curl(vectorfield):

A=Magnetic_Potenial_field(2,4,3,0, 1,.5,1 ,10, 10,10)
for i in range(np.shape(A[0])[1]):
    for j in range(np.shape(A[0])[0]):
        print("x:",A[0][i][j],"y:",A[1][i][j],"z:",A[2][i][j],end=" ")
    print()
