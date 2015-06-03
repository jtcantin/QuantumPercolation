import sys
import numpy as np
import matplotlib.pyplot as plt

eigvalFilename = sys.argv[1]
nu = int(sys.argv[2])
nBins = int(sys.argv[3])

eigvals = np.loadtxt(eigvalFilename,delimiter=",")

#print eigvals.size
#print eigvals

#Calculate the local level density, as a function of nu
locLevelRhoLarge = -2*nu / (np.roll(eigvals,nu) - np.roll(eigvals,-nu))

#Ignore the end values as they don't have an appropriate value
locLevelRho = locLevelRhoLarge[nu:(-nu)]

#Get the energy level spacings
spacings = np.diff(eigvals)

#Normalize the spacing with the local level density
spacingsNormalized = (spacings*locLevelRhoLarge[:-1])[nu:(-nu)]

#Bin the data
binnedData, bin_edges = np.histogram(spacingsNormalized, bins=nBins)

#print bin_edges
#print binnedData

#Calculate Level density from 2nd order derivative
#O2rho = 1. / np.gradient(eigvals


#Plot
figNum = 0

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
line, = ax.plot(eigvals[nu:(-nu)], locLevelRho)

#plt.ylim(0,0.1)
#plt.xlim(0,N)
plt.xlabel("Energy", fontsize=16)
plt.ylabel(r'$\rho$', fontsize=16)
plt.title("Local Level Density", fontsize=18)

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
bar = ax.bar(bin_edges[:-1], binnedData,width=bin_edges[1]-bin_edges[0])

#plt.ylim(0,0.1)
#plt.xlim(0,N)
plt.xlabel(r'$s$', fontsize=16)
plt.ylabel(r'$P(s)$', fontsize=16)
plt.title("Energy Level Statistics", fontsize=18)

plt.show()



















