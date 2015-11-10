#To run: python LDOS_IPRCalc.py [evecFilename] [evalFilename] [evecFilename] [evalFilename] [evecFilename] [evalFilename] ...

import sys
import numpy as np
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

np.set_printoptions(threshold=10000,linewidth=2000,precision=3,suppress=False)

Ndis = (len(sys.argv) - 1)/2.
print "Number of Disorders:", Ndis

#binnedData_Avg = np.zeros(1)
#bin_edges_Avg = ??

E_Array = np.linspace(-10.0, 10.0, 1000)
LDOS_Avg_Array_Avg = np.zeros(E_Array.size)
LDOS_Typ_Array_Avg_PreExp = np.zeros(E_Array.size)

def gauss(x,mu,sigma,denom):
    exponent = -1. * ((x - mu)**2) / (2.*(sigma**2))
    return denom * np.exp(exponent)

def gaussArray(x,mu,sigma,denom,sqrLapackEvec_Array):
    gaussVal = np.zeros(x.size)
    for i in range(0,mu.size):
        exponent = -1. * ((x - mu[i])**2) / (2.*(sigma**2))
        gaussVal += sqrLapackEvec_Array[i] * denom * np.exp(exponent)
    return gaussVal

#sigma = (E_Array[1] - E_Array[0])/2.
sigma = 500.
denom = (1./(sigma*np.sqrt(2.*np.pi)))

for disOrderNum in range(1, len(sys.argv), 2):
    LapackEvecFilename = sys.argv[disOrderNum]
    LapackEvalFilename = sys.argv[disOrderNum + 1]


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
    LDOS_Broadened = gaussArray(E_Array,Eval_array,sigma,denom,sqrLapackEvec_Array)
    LDOS_Avg_Array = np.sum(LDOS_Broadened, axis = 0) / basisSize
    LDOS_Typ_Array = np.sum(np.log(LDOS_Broadened), axis=0) / basisSize



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

    #Average over values
#    IPR_Array_List.append(
#    binnedData_Avg += binnedData / Ndis
#    bin_edges_Avg = ??
    LDOS_Avg_Array_Avg += LDOS_Avg_Array / Ndis
    LDOS_Typ_Array_Avg_PreExp += LDOS_Typ_Array / Ndis
#    Eval_array_Avg += Eval_array / NDis

LDOS_Typ_Array_Avg = np.exp(LDOS_Typ_Array_Avg_PreExp)

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

#line = ax.plot(E_Array, LDOS_Avg_Array_Avg, 'r-', label="Average LDOS")
line2 = ax.plot(E_Array, LDOS_Typ_Array_Avg, 'k-', label="Typical LDOS")
line3 = ax.plot(E_Array, LDOS_Broadened, 'b-', label="Broadened LDOS")

ax.legend()

plt.xlabel("DOS", fontsize=16)
plt.ylabel("Energy", fontsize=16)



plt.show()







