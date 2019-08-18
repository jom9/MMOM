import numpy as np
class Optical_Functions():

    def __init__(self,minWavelength,maxWaveLength,step):
        self.num_of_points=0
        self.minWavelength=minWavelength
        self.maxWaveLength=maxWaveLength
        lamda_incremeter=minWavelength
        self.wavelength=[]
        while lamda_incremeter<maxWaveLength:
            self.wavelength.append(lamda_incremeter)
            lamda_incremeter+=step
class SiO2(Optical_Functions):

    def get_n(self):
        n=[]
        for lamda in self.wavelength:
            if lamda>1.5 and lamda<14.5:
                n.append(self.Kischkat_n(lamda))
            elif lamda<=1.5:
                n.append(self.Rodriguez_n(lamda))
            else:
                raise Outside_Wavelength_Range
        return n
    def get_k(self):
        k=[]
        for lamda in self.wavelength:
            if lamda>1.5 and lamda<14.5:
                k.append(self.Kischkat_k(lamda))
            elif lamda<=1.5:
                k.append(self.Rodriguez_k(lamda))
            else:
                raise Outside_Wavelength_Range
        return k
    def Kischkat_n(self,lamda):
        #https://refractiveindex.info/?shelf=main&book=SiO2&page=Kischkat
        Numcoeffs =   [1.529,-56.89,821.2,-5661,1.806*10**4,-1.975*10**4]
        Demcoeffs =[1,-37.62,549,-3823,1.231*10**4,-1.355*10**4]
        if lamda>1.5 and lamda<14.5:
                i=0
                n=0

                for coeff in Numcoeffs:
                    n+=coeff*(lamda)**(len(Numcoeffs)-i)
                    i+=1
                d=0
                i=0
                for coeff in Demcoeffs:
                    d+=coeff*(lamda)**(len(Demcoeffs)-i)
                    i+=1
                return n/d
        else:
            raise Outside_Wavelength_Range
    def Kischkat_k(self,lamda):
        #https://refractiveindex.info/?shelf=main&book=SiO2&page=Kischkat
        Numcoeffs =   [-0.03271,0.6908,-3.989,8.193,-5.497]
        Demcoeffs = [1,-20.4,121.2,-170,43.37]
        if lamda>1.5 and lamda<14.5:
                i=0
                k=0

                for coeff in Numcoeffs:
                    k+=coeff*(lamda)**(len(Numcoeffs)-i)
                    i+=1
                d=0
                i=0
                for coeff in Demcoeffs:
                    d+=coeff*(lamda)**(len(Demcoeffs)-i)
                    i+=1
                return k/d
        else:
            raise Outside_Wavelength_Range
    def Rodriguez_n(self,lamda):
        #https://refractiveindex.info/?shelf=main&book=SiO2&page=Rodriguez-de_Marcos
        Numcoeffs =   [4068,8308,9246,5104,-1195,76.08]
        Demcoeffs =[1,1.47*10**4,-3345 ,7706,-1578,93.91]



        if lamda>0 and lamda<1.5:
                i=0
                n=0

                for coeff in Numcoeffs:
                    n+=coeff*(lamda)**(len(Numcoeffs)-i)
                    i+=1
                d=0
                i=0
                for coeff in Demcoeffs:
                    d+=coeff*(lamda)**(len(Demcoeffs)-i)
                    i+=1
                return n/d
        else:
            raise Outside_Wavelength_Range
    def Rodriguez_k(self,lamda):
        #https://refractiveindex.info/?shelf=main&book=SiO2&page=Rodriguez-de_Marcos
        Numcoeffs =   [-0.009503,0.00895,-0.002131,0.0001589,-1.257*10**-6]
        Demcoeffs = [1,-0.3789, 0.05341,-0.003366,8.336*10**-5]

        if lamda>0 and lamda<1.5:
                i=0
                k=0

                for coeff in Numcoeffs:
                    k+=coeff*(lamda)**(len(Numcoeffs)-i)
                    i+=1
                d=0
                i=0
                for coeff in Demcoeffs:
                    d+=coeff*(lamda)**(len(Demcoeffs)-i)
                    i+=1
                return k/d
        else:
            raise Outside_Wavelength_Range

class Air(Optical_Functions):
    def get_n(self):
        n=[]
        for lamda in self.wavelength:
            n+=[1]
        return n
    def get_k(self):
        k=[]
        for lamda in self.wavelength:
            k+=[0]
        return k
class PMMA(Optical_Functions):
    def get_n(self):
        n=[]

        for lamda in self.wavelength:
            if lamda<2.5:
                n.append(self.Tsuda_n(lamda))
            elif lamda>=2.5 and lamda<6:
                n.append(self.Tsuda_Drude_n1(lamda))
            elif lamda>=6 and lamda<18:
                n.append(self.Tsuda_Drude_n2(lamda))
            else:
                raise Outside_Wavelength_Range
        return n

    def get_k(self):
        k=[]
        for lamda in self.wavelength:
            if lamda<2.5:
                k.append(self.Tsuda_k(lamda))
            else:
                k.append(self.Tsuda_Drude_k(lamda))
        return k
    def Tsuda_n(self,lamda):
        #https://refractiveindex.info/?shelf=organic&book=poly%28methyl_methacrylate%29&page=Tsuda
        if lamda<2.5 and lamda>=0:
            return 1.470+.008354/(lamda**2)-.0008309/(lamda**4)
        else:
            raise Outside_Wavelength_Range
    def Tsuda_k(self,lamda):
        #https://refractiveindex.info/?shelf=organic&book=poly%28methyl_methacrylate%29&page=Tsuda
        if lamda<2.5 and lamda>=0:
            return 0
        else:
            raise Outside_Wavelength_Range
    def Tsuda_Drude_n1(self,lamda):
        if lamda<=0 and lamda<6:
            Numcoeffs =   [ 1.345,-13.87,26.08,52.17,2.147,19.33]
            Demcoeffs =[1,-11.05,28.98 ,0.3563,47.56,-1.696]
            i=0
            n=0
            for coeff in Numcoeffs:
                n+=coeff*(lamda)**(len(Numcoeffs)-i)
                i+=1
            d=0
            i=0
            for coeff in Demcoeffs:
                d+=coeff*(lamda)**(len(Demcoeffs)-i)
                i+=1
            return n/d
        else:
            raise Outside_Wavelength_Range
    def Tsuda_Drude_n2(self,lamda):
        if lamda>6 and lamda<18:
            n=1.51
            Acoeffs = [ 0.03275 ,-0.01605,-0.01463,0.01794,-0.002162,-0.01627,0.01049,0.01528]
            Bcoeffs = [0.006733,0.004616,-0.02704,0.02429,-0.0009015,-0.01631,0.001763, 0.01916,-0.01087 ]
            w =0.5174
            for i in range(len(Acoeffs)):
                n+=np.cos(w*Acoeffs[i])+np.sin(w*Bcoeffs[i])
            return n
        else:
            raise Outside_Wavelength_Range
    def Tsuda_Drude_k(self,lamda):
        Acoeffs= [0.8239,0.3107,0.2039,0.1197,0.072, 0.03177, 0.04899]
        Bcoeffs= [5.775,8.686, 8.396,8.005,6.884,9.486, 13.3]
        Ccoeffs= [0.02582,0.1562,0.09626,0.2097, 0.1536 ,3.094 ,0.1523]
        if lamda>2.5 and lamda<18:
            n=0
            for i in range(len(Acoeffs)):
                n+=Acoeffs[i]*np.exp(   -((lamda-Bcoeffs[i])/Ccoeffs[i])**2  )
            return n
        else:
            raise Outside_Wavelength_Range
class PVA(Optical_Functions):
    def get_n(self):
        n=[]
        for lamda in self.wavelength:
            n+=[1.460+0.00665/(lamda**2)]
        return n
    def get_k(self):
        k=[]
        for lamda in self.wavelength:
            k+=[0]
        return k
class F8BT(Optical_Functions):
    def get_n(self):
        nth=[]

        for lamda in self.wavelength:
            nth.append(self.n(lamda))
        return nth
    def get_k(self):
        k=[]
        for lamda in self.wavelength:
            k.append(self.k(lamda))
        return k
    def n(self,lamda):
        n=0
        Acoeffs= [0.4828,0.2266,0.3468,1.81]
        Bcoeffs= [0.5029,0.5953,0.3512,1.29]
        Ccoeffs=[0.04006,0.1049, 0.02751,2.007]
        for i in range(len(Acoeffs)):
            n+=Acoeffs[i]*np.exp(   -((lamda-Bcoeffs[i])/Ccoeffs[i])**2  )
        return n
    def k(self,lamda):
        k=0
        Acoeffs= [0,0,0.3058,0.415,0.4212,0.1292]
        Bcoeffs= [ 0.4749,0.4743,0.4801,0.443, 0.3249 ,1.076]
        Ccoeffs=[0.0004661,2.22*10**-14,0.03027,0.0444,0.03442,0.4018]
        for i in range(len(Acoeffs)):
            k+=Acoeffs[i]*np.exp(   -((lamda-Bcoeffs[i])/Ccoeffs[i])**2  )
        return k
class CA(Optical_Functions):
    def get_n(self):
        n=[]
        for lamda in self.wavelength:
            n+=[self.n(lamda)]
        return n
    def get_k(self):
        n=[]
        for lamda in self.wavelength:
            n+=[self.k(lamda)]
        return n
    def n(self,lamda):
            return 1.534+0.03305/(lamda**2)
    def k(self,lamda):
        n=2.1979*10**(-5)
        Acoeffs = [6.817*10**-6,2.626*10**-6]
        Bcoeffs = [1.235*10**-5,-5.54*10**-6]
        w =       6.787
        for i in range(len(Acoeffs)):
            n+=np.cos(w*Acoeffs[i]*lamda)+np.sin(w*Bcoeffs[i]*lamda)
        return n
class Outside_Wavelength_Range(Exception):
    "No fitting function in that range"
    pass
