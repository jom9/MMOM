from Magnet import Magnet
import numpy as np
import os
class MARS():

    def __init__(self,xlim,ylim,zlim,resolution,timesres,timelim,step):
        self.xlim=xlim
        self.ylim=ylim
        self.zlim=zlim
        self.resolution=resolution
        self.timesres=timesres
        self.step=step

class Simple_Gear(MARS):
    def set_parameters(self,inner_num_of_teeth,teeth_ratio,thickness,diameter,inner_armLength,outer_armlength,magnetization,driving_angular_velocity):
        self.inner_num_of_teeth=inner_num_of_teeth
        self.teeth_ratio=teeth_ratio
        self.thickness=thickness
        self.diameter=diameter
        self.inner_armLength=inner_armLength
        self.outer_armlength=outer_armlength
        self.magnetization=magnetization
        self.driving_angular_velocity=driving_angular_velocity
    def genDrivingBfield(self):
        Driver =Magnet(self.thickness,self.diameter,self.inner_armLength,self.magnetization,self.xlim,self.ylim,self.zlim,self.resolution,self.step)
        phi=0

        for i in range(int((2*np.pi/self.driving_angular_velocity)/self.timesres) ):
            Driver.genBfield(phi,i)
            phi= phi+i*self.timesres*self.driving_angular_velocity


        for i in range(int(   (2*np.pi/self.driving_angular_velocity)/self.timesres)):
            Bfield_list=[]


            for k in range(self.inner_num_of_teeth):
                j=int((int(   (2*np.pi/self.driving_angular_velocity)/self.timesres)/self.inner_num_of_teeth)*k)
                if i+j<int((2*np.pi/self.driving_angular_velocity)/self.timesres):

                    Bfield_list.append(Driver.getBfield(i+j))
                else:
                    reducediplusj=i+j
                    while reducediplusj>=int((2*np.pi/self.driving_angular_velocity)/self.timesres):
                        reducediplusj-=int((2*np.pi/self.driving_angular_velocity)/self.timesres)
                    Bfield_list.append(Driver.getBfield(reducediplusj))
            Bx_field = np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
            By_field = np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
            Bz_field = np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
            while Bfield_list:
                B_field=Bfield_list.pop()
                Bx_field=np.add(Bx_field,B_field[0])
                By_field=np.add(Bx_field,B_field[1])
                Bz_field=np.add(Bx_field,B_field[2])
            current_path = os.getcwd()
            if "driver_Bfielddata" not in os.listdir(current_path):
                os.mkdir("driver_Bfielddata")
            os.chdir("driver_Bfielddata")
            file= open("driver_Bfielddata"+str(i)+".txt",'w+')

            for u in range(int(2*self.xlim/self.resolution)):
                for v in range(int(2*self.ylim/self.resolution)):
                    for w in range(int(2*self.zlim/self.resolution)):
                        file.write( "("+str(u*self.resolution-self.xlim)+","+str(v*self.resolution-self.ylim)+","+str(w*self.resolution-self.zlim)+")  "+"("+str(Bx_field[u][v][w])+","+str(By_field[u][v][w])+","+str(Bz_field[u][v][w]) +")\r\n")
            os.chdir(current_path)

    def get_Driving_B_field(self,i):
        try:
            current_path=os.getcwd()
            os.chdir("driver_Bfielddata")
            file = open("driver_Bfielddata"+str(i)+".txt",'r')
            Bx_field = np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
            By_field = np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
            Bz_field = np.zeros((int(2*self.xlim/self.resolution),int(2*self.ylim/self.resolution),int(2*self.zlim/self.resolution)))
            for line in file.readlines():
                parsedline = line.split(" ")

                coors = parsedline[0].split(",")
                field = parsedline[-1].split(",")
                if len(coors)==3 and len(field)==3:


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
