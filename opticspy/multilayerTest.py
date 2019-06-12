from multilayers import multilayers
from importer import importer

f=importer()

(n,d)=f.importMultiLayers('multilayertestdata.txt')
lamda = 5*10**-6
print( multilayers(n,d).Transmittance(lamda))
