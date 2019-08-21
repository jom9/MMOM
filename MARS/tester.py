import MARS
import numpy as np
xlim=15
ylim=15
zlim=15
res= .5
timesres=1
timelim=72
step=.5
Gear = MARS.Simple_Gear(xlim,ylim,zlim,res,timesres,timelim,step)

inner_num_of_teeth=4
teeth_ratio=4

thickness=4
diameter=2
inner_armLength=5

outer_armlength=10

magnetization=1
angular_velocity=2*np.pi/72

Gear.set_parameters(inner_num_of_teeth,teeth_ratio,thickness,diameter,inner_armLength,outer_armlength,magnetization,angular_velocity)
#Gear.genDrivingBfield()
