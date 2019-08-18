from waveguides import waveguide
import functions
import matplotlib.pyplot as plt
import numpy as np
lamda1=.850
lamda2=1.300
lamda3=1.550

thickness = []
width =10**-5


amplitude=1
c=3*10**8
thickness= np.arange(.5,1.5,.01)
loss1=0*thickness
loss2=0*thickness
loss3=0*thickness
phase=0
i=0
for t in thickness:
    frequency=2*np.pi*c/(lamda1*(10**-9))
    ntop=functions.PVA(.1,2,.01).get_n()[84]-.2
    ncenter=functions.Air(.1,2,.01).get_n()[84]
    nbottom=functions.PVA(.10,2,.01).get_n()[84]-.2
    ktop=functions.PVA(.10,2,.01).get_k()[84]
    kcenter=functions.Air(.1,2,.01).get_k()[84]
    kbottom=functions.PVA(.10,2,.01).get_k()[84]
    w=waveguide(amplitude,phase,ntop,ncenter,nbottom,ktop,kcenter,kbottom,frequency,t,width)
    loss1[i]=1.25*w.loss('s',.00025*np.pi)
    print(t,ntop,ktop,loss1[i])
    frequency=2*np.pi*c/(lamda2*(10**-9))
    ntop=functions.PVA(.10,2,.01).get_n()[129]-.2
    ncenter=functions.Air(.1,2,.01).get_n()[129]
    nbottom=functions.PVA(.10,2,.01).get_n()[129]-.2
    ktop=functions.PVA(.10,2,.01).get_k()[129]
    kcenter=functions.Air(.1,2,.01).get_k()[129]
    kbottom=functions.PVA(0.1,2,.01).get_k()[129]
    w=waveguide(amplitude,phase,ntop,ncenter,nbottom,ktop,kcenter,kbottom,frequency,t,width)
    loss2[i]=1.1*w.loss('s',.00025*np.pi)
    print(t,loss2[i])
    frequency=2*np.pi*c/(lamda3*(10**-9))
    ntop=functions.PVA(0.1,2,.01).get_n()[154]-.2
    ncenter=functions.Air(.1,2,.01).get_n()[154]
    nbottom=functions.PVA(0.1,2,.01).get_n()[154]-.2
    ktop=functions.PVA(0.1,2,.01).get_k()[154]
    kcenter=functions.Air(0.1,2,.01).get_k()[154]
    kbottom=functions.PVA(0.1,2,.01).get_k()[154]
    w=waveguide(amplitude,phase,ntop,ncenter,nbottom,ktop,kcenter,kbottom,frequency,t,width)
    loss3[i]=w.loss('s',.00025*np.pi)
    print(t,loss3[i])
    i+=1

plt.subplot(111)
plt.plot(thickness,loss1)
plt.title("Signal Loss in PANI Waveguide with Air Core at 850 nm")
plt.xlabel("Length (cm)")
plt.ylabel("Loss dB")
plt.show()

plt.subplot(111)
plt.plot(thickness*100,loss2)
plt.title("Signal Loss in PANI Waveguide with Air Core at 1300 nm")
plt.xlabel("Length (cm)")
plt.ylabel("Loss dB")
plt.show()

plt.subplot(111)
plt.plot(thickness*100,loss3)
plt.title("Signal Loss in PANI Waveguide with Air Core at 1550 nm")
plt.xlabel("Length (cm)")
plt.ylabel("Loss dB")

plt.show()

plt.subplot(111)
plt.plot(thickness*100,loss1)
plt.plot(thickness*100,loss2)
plt.plot(thickness*100,loss3)
plt.title("Signal Loss in PANI Waveguide with Air Core")
plt.xlabel("Length (cm)")
plt.ylabel("Loss dB")
plt.legend(['850nm', '1300nm', '1500nm'])

plt.show()
