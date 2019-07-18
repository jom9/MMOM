from sympy import *
import numpy as np
import time
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

    cornerPos= currentPos(thickness,diameter,armLength,phi) #gets the four corners needed to take each rectangluar slice
    posFunction = postionFunctions(thickness,diameter,armLength,phi) #creates the lines that trace out borders to slice

    Ax = numeric_double_Integral(integrand[0],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Ax += numeric_double_Integral(integrand[0],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Ax += numeric_double_Integral(integrand[0],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)

    Ay = numeric_double_Integral(integrand[1],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Ay += numeric_double_Integral(integrand[1],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Ay += numeric_double_Integral(integrand[1],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)

    Az = numeric_double_Integral(integrand[2],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Az += numeric_double_Integral(integrand[2],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
    Az += numeric_double_Integral(integrand[2],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)
    # evualtes the integral to get the potential. each component takes 3 integrals since the region must be broken up into 3 subregions
    A=np.array([Ax,Ay,Az])

    return A

def Magnetic_Potenial_field(thickness,diameter,armLength,phi, magnetization,step,resolution ,xlim, ylim,zlim):
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
        k+=step
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
                print("stuff",i,j,k)
                AijkX[i][j][k]=A(i*resolution-xlim,j*resolution-ylim,k*resolution-zlim)[0]
                AijkY[i][j][k]=A(i*resolution-xlim,j*resolution-ylim,k*resolution-zlim)[1]
                AijkZ[i][j][k]=A(i*resolution-xlim,j*resolution-ylim,k*resolution-zlim)[2]
                k+=1
            j+=1
        i+=1
    return np.array([AijkX,AijkY,AijkZ])
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

    return sortbyx([r0,r1,r2,r3])
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
    print((round(L[0][0],2),round(L[0][1],2) ),(round(L[1][0],2),round(L[1][1],2) ),(round(L[2][0],2),round(L[2][1],2) ),(round(L[3][0],2),round(L[3][1],2) ) )
    print("current pos",(round(slopeofLine0To1,2),round(b01,2)),(round(slopeofLine0To2,2),round(b02,2)),(round(slopeofLine1To3,2),round(b13,2)),(round(slopeofLine2To3,2),round(b23,2)))
    return [(slopeofLine0To1,b01),(slopeofLine0To2,b02),(slopeofLine1To3,b13),(slopeofLine2To3,b23)]
def numeric_double_Integral(integrand,x0,xf,yb,yt,xstep,ystep,symx,symy):
    #simple rieman double integral
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
        numberofticks= int(totaltime/timestep)
        Bfield=[] # bfield
        zprime =0
        for t in range(numberofticks):
            print("started",t,"!")
            start = time.time()
            phi=angularvelocity*t+startphi

            M = np.array([np.cos(phi),np.sin(phi),0]) # magnetization vector

            x,y,z,xprime,yprime= symbols('x y z xprime yprime')
            delR = np.array([x-xprime,y-yprime,z-zprime]) # this is the distance vector, with the primed values denoting source points and the unprimed values denoting field values

            integrand =['1','1','1']#np.cross(M,delR)/(magnitude(delR)**3)

            cornerPos= currentPos(thickness,diameter,armLength,phi) #gets the four corners needed to take each rectangluar slice
            posFunction = postionFunctions(thickness,diameter,armLength,phi) #creates the lines that trace out borders to slice

            Ax = numeric_double_Integral(integrand[0],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            Ax += numeric_double_Integral(integrand[0],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            Ax += numeric_double_Integral(integrand[0],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)

            Ay = numeric_double_Integral(integrand[1],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            Ay += numeric_double_Integral(integrand[1],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            Ay += numeric_double_Integral(integrand[1],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)

            Az = numeric_double_Integral(integrand[2],cornerPos[0][0],cornerPos[1][0],xprime*posFunction[0][0]+posFunction[0][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            Az += numeric_double_Integral(integrand[2],cornerPos[1][0],cornerPos[2][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[1][0]+posFunction[1][1],step,xprime,yprime)
            Az += numeric_double_Integral(integrand[2],cornerPos[2][0],cornerPos[3][0],xprime*posFunction[2][0]+posFunction[2][1],xprime*posFunction[3][0]+posFunction[3][1],step,xprime,yprime)
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
            A=Magnetic_Potenial_field(thickness,diameter,armLength,phi, magnetization,step,resolution ,xlim, ylim,zlim)
            #Bfield.append(Magnetic_Field(A,step))
            file = open("Bfielddata"+str(t)+".txt","w+")
            B=Magnetic_Field(A,res)
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

res =1
xlim =10
ylim =10
zlim =10
step=.5
totaltime=50
timestep=1
startphi =0*np.pi/2
thickness=4
diameter=2
armLength=3
angularvelocity=np.pi/3
magnetization=200
#rotating_field(thickness,diameter,armLength,startphi, angularvelocity,magnetization,step,resolution ,xlim, ylim,zlim, totaltime,timestep):
double_integral_checker(thickness,diameter,armLength,startphi, angularvelocity,magnetization,step,res ,xlim, ylim,zlim, totaltime,timestep)
#B =rotating_field(thickness,diameter,armLength,startphi, angularvelocity,magnetization,step,res ,xlim, ylim,zlim, totaltime,timestep)
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
