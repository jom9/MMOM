import numpy as np
import tkinter

class importer:

	def __init__(self):
		pass
	def importIndexOfRefraction(self,filename):
		file = open(filename,'r')
		lines = file.readlines()
		numofdatapoints = len(lines)
		n = np.zeros((numofdatapoints))
		lamnda = np.zeros((numofdatapoints))
		i=0
		for line in lines:
			lamda[i]= line.split()[0]
			n[i] = line.split()[1]
			i+=1
		return (lamda,n)

	def importComplexIndexOfRefraction():
		file = open(filename,'r')
		lines = file.readlines()
		numofdatapoints = len(lines)
		k = np.zeros((numofdatapoints))
		lamnda = np.zeros((numofdatapoints))
		i=0
		for line in lines:
			lamda[i]= line.split()[0]
			k[i] = line.split()[1]
			i+=1
		return (lamda,k)

	def importMultiLayers(self,filename):
		file = open(filename,'r')
		lines = file.readlines()
		numofdatapoints = len(lines)
		d = np.zeros((numofdatapoints)) * 10**-6 #thickness is the micrometer range
		n = np.zeros((numofdatapoints))
		i=0
		for line in lines:
			d[i]= line.split()[0]
			n[i] = line.split()[1]
			i+=1
		return (n,d)
