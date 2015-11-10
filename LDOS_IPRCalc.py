#To run: python LDOS_IPRCalc.py evecFilename evalFilename

import sys
import numpy as np
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

np.set_printoptions(threshold=10000,linewidth=2000,precision=3,suppress=False)

LapackEvecFilename = sys.argv[1]
LapackEvalFilename = sys.argv[2]


LapackEvecFile = open(LapackEvecFilename, 'r')

LapackEvec_List = []

#Skip first two comment lines
LapackEvecFile.readline()
LapackEvecFile.readline()

for line in LapackEvecFile:
    splitLine = line.split()
    strip_list = map(lambda it: it.strip(), splitLine)
    LapackEvec_List.append(np.array(map(float, strip_list[1:])))

LapackEvecFile.close()

LapackEvec_Array = np.vstack(LapackEvec_List)
#print LapackEvec_Array
#print LapackEvec_Array[:,0]
#print LapackEvec_Array[0]
#print LapackEvec_Array.shape
basisSize = LapackEvec_Array[:,0].size
numEigStates = LapackEvec_Array[0].size

print "------------------------"
LapackEvecFilename_short = LapackEvecFilename.split("/")[-1]
print "For file ", LapackEvecFilename_short, ": "
print "Basis Size: ", basisSize
print "Number of Eigenstates in file: ", numEigStates

#Check orthonormality
#print np.dot(LapackEvec_Array.T,LapackEvec_Array)
absDiffOrtho = np.abs(np.identity(numEigStates) - np.dot(LapackEvec_Array.T,LapackEvec_Array))
#print absDiffOrtho
maxAbsDiff = np.max(absDiffOrtho)
print "Max Abs Diff from Orthonormality: ", maxAbsDiff

#Calculate the IPR I_2,nu as define in eqn. 2 of DOI: 10.1103/PhysRevB.83.184206
#Square the coefficients:
sqrLapackEvec_Array = LapackEvec_Array**2
numerator_Array = np.sum(sqrLapackEvec_Array**2, axis=0)
denomentator_Array = np.sum(sqrLapackEvec_Array,axis=0)**2
IPR_Array = numerator_Array/denomentator_Array

#Calculate Average and Typical DOS

#Import Eigenvalues
Eval_array = np.loadtxt(LapackEvalFilename,delimiter=",")
#print Eval_array

#LapackEvalFile = open(LapackEvalFilename, 'r')
#
#LapackEval_List = []
#
#for line in LapackEvalFile:
#    LapackEval_List.append(float(line))
#
#LapackEvalFile.close()

#Calculate ADOS and TDOS
LDOS_Avg_Array = np.sum(sqrLapackEvec_Array, axis = 0) / basisSize
LDOS_Typ_Array = np.exp(np.sum(np.log(sqrLapackEvec_Array), axis=0) / basisSize)



#print LapackEvec_Array
#print sqrLapackEvec_Array
#print numerator_Array
#print denomentator_Array
#print IPR_Array
print "Min IPR: ", np.min(IPR_Array)
print "Max IPR: ", np.max(IPR_Array)

#Arrange into a histogram
nBins = 100
binnedData, bin_edges = np.histogram(IPR_Array, bins=nBins,density=True)

#########################
#########################
#Plot
#########################
#########################

figNum = 0

#IPR Histogram
figNum += 1
fig = plt.figure(figNum,facecolor="white")
ax = plt.subplot()
bar = ax.bar(bin_edges[:-1], binnedData,width=bin_edges[1]-bin_edges[0])
#bar = ax.bar(bin_edges[1:-1], binnedData[1:],width=bin_edges[1]-bin_edges[0])

#plt.ylim(0,0.1)
plt.xlim(0,1)

#Add vertical line at 1/N, where N is number of sites
print "Ideal Fully Delocalized IPR: ", 1./basisSize
#print ax.get_ylim()
#line = plt.plot([1./basisSize, 1./basisSize], [0,ax.get_ylim()[1]], 'r', linewidth=2)


#line = plt.plot([0.005, 0.005], [0,3000], 'k')

plt.xlabel("IPR", fontsize=16)
plt.ylabel("Probability Density", fontsize=16)
#plt.title("IPR Histogram", fontsize=18)


filenameList = LapackEvecFilename.split("/")
FigureName = filenameList[-2] + "_"+ filenameList[-1][13:-12] + "_IPRdist.png"
#FigureName = filenameList[-1][:-4] + "_IPRdist.png"
#FigureNameFile = "/".join(filenameList[:-1]+[FigureName])
#FigureNameFile = "/home/jtcantin/code/QuantumPercolation/srcTRX/data/IPRGraphs/" + FigureName
FigureNameFile = "/Volumes/Fano/code/QuantumPercolation/srcTRX/data/IPRGraphs/" + FigureName


#print filenameList
#print FigureName
#print FigureNameFile

#plt.savefig(FigureNameFile, bbox_inches='tight')

#IPR as a function of energy
figNum += 1
fig = plt.figure(figNum,facecolor="white")
ax = plt.subplot()

line = ax.plot(Eval_array, IPR_Array, 'k-')

#ax.legend()

plt.ylabel("IPR", fontsize=16)
plt.xlabel("Energy", fontsize=16)



#LDOS Plot as a function of energy
figNum += 1
fig = plt.figure(figNum,facecolor="white")
ax = plt.subplot()

#line = ax.plot(Eval_array, LDOS_Avg_Array, 'r-', label="Average LDOS")
line2 = ax.plot(Eval_array, LDOS_Typ_Array, 'k-', label="Typical LDOS")

ax.legend()

plt.ylabel("DOS", fontsize=16)
plt.xlabel("Energy", fontsize=16)



plt.show()







