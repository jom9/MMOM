from Magnet import Magnet
import numpy as np
class MARS():

    def __init__(self,xlim,ylim,zlim,resolution,timesres,timelim,):
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
    def gen_Driving_B_field(self):
        Driver =Magnet(self.thickness,self.diameter,self.inner_armLength,self.xlim,self.ylim,self.zlim,self.resolution,self.step)
        phi=0
        for i in range(int((2*np.pi/self.driving_angular_velocity)/self.timesres) ):
            Driver.gen_B_field(phi,i)
            phi= phi+i*self.timesres*self.driving_angular_velocity

        for i in range(int((2*np.pi/self.driving_angular_velocity)/self.timesres)):
            Bfield_list=[]
            file= open("driver_Bfielddata"+str(i)+".txt",'w+')

            for j in range(inner_num_of_teeth):
                if i+j<int((2*np.pi/self.driving_angular_velocity)/self.timesres):
                    Bfield_list.append(Driver.get_B_field(i+j))
                else:
                    reducediplusj=i+j
                    while i+j>=int((2*np.pi/self.driving_angular_velocity)/self.timesres):
                        reducediplusj-=int((2*np.pi/self.driving_angular_velocity)/self.timesres)
                    Bfield_list.append(Driver.get_B_field(reducediplusj))
            Bx_field = np.zeros((2*self.xlim/self.resolution,2*self.ylim/self.resolution,2*self.zlim/self.resolution))
            By_field = np.zeros((2*self.xlim/self.resolution,2*self.ylim/self.resolution,2*self.zlim/self.resolution))
            Bz_field = np.zeros((2*self.xlim/self.resolution,2*self.ylim/self.resolution,2*self.zlim/self.resolution))
            while Bfield_list:
                B_field=Bfield_list.pop()
                Bx_field=np.add(Bx_field,B_field[0])
                By_field=np.add(Bx_field,B_field[1])
                Bz_field=np.add(Bx_field,B_field[2])
            for u in range(2*self.xlim/self.resolution):
                for v in range(2*self.ylim/self.resolution):
                    for w in range(2*self.zlim/self.resolution):
                        file.write( "("+str(u*self.resolution-self.xlim)+","+str(v*self.resolution-self.ylim)+","+str(w*self.resolution-self.zlim)+")  "+"("+str(Bx_field[u][v][w])+","+str(By_field[u][v][w])+","+str(Bz_field[u][v][w]) +")\r\n")


    def get_Driving_B_field(self,i):
        file = open("driver_Bfielddata"+str(i)+".txt",'r')
        Bx_field = np.zeros((2*self.xlim/self.resolution,2*self.ylim/self.resolution,2*self.zlim/self.resolution))
        By_field = np.zeros((2*self.xlim/self.resolution,2*self.ylim/self.resolution,2*self.zlim/self.resolution))
        Bz_field = np.zeros((2*self.xlim/self.resolution,2*self.ylim/self.resolution,2*self.zlim/self.resolution))
        for line in file.readlines():
            parsedline = line.split(" ")
