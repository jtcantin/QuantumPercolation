# To run: python eLevelStats.py [eigvalFilename] [nu] [nBins] [sigma]

import sys
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(threshold=10000,linewidth=2000,precision=20,suppress=False)

eigvalFilename = sys.argv[1]
nu = int(sys.argv[2])
nBins = int(sys.argv[3])
sigma = float(sys.argv[4])

eigvals = np.loadtxt(eigvalFilename,delimiter=",")
N = eigvals.size
#print eigvals.size
#print eigvals[:10]

#Calculate the local level density, as a function of nu
locLevelRhoLarge = -2*nu / (np.roll(eigvals,nu) - np.roll(eigvals,-nu))

#Ignore the end values as they don't have an appropriate value
locLevelRho = locLevelRhoLarge[nu:(-nu)]

print "Local level density calculated."

#Get the energy level spacings
spacings = np.diff(eigvals)

#Normalize the spacing with the local level density
spacingsNormalized = (spacings*locLevelRhoLarge[:-1])[nu:(-nu)]

#Bin the data
binnedData, bin_edges = np.histogram(spacingsNormalized, bins=nBins)

print "Local level density binned."

#########################
#Calculate the consistency requirement (see eqn 4.8.2 of Quantum Signatures of Chaos, 2010 by Haake):
#########################

locLevelRhoDer = (np.roll(locLevelRho,-1)[1:-1] - np.roll(locLevelRho,1)[1:-1]) / ( np.abs(np.diff(eigvals[nu:(-nu)])) + np.abs(np.roll(np.diff(eigvals[nu:(-nu)]),1)) )[1:]

#print locLevelRhoDer
checkRatio = np.abs(locLevelRhoDer) / (locLevelRho[1:-1]**2)
#print checkRatio
print "Maximum Consistency Ratio (want << 1): ", np.nanmax(checkRatio)


#print bin_edges
#print binnedData
#print np.mean(binnedData)
#print np.max(binnedData)
#print np.min(binnedData)
#print binnedData[np.nonzero(binnedData)]

#########################
#Try an alternative method: NOT WORKING
#   Integrate over rho and take the differences
#########################
#locLevelRhoInt = locLevelRho*(((np.roll(eigvals,-1)-np.roll(eigvals,1))[nu:(-nu)])/2.)
#
#NEspacings = np.diff(locLevelRhoInt)
#binnedData2, bin_edges2 = np.histogram(NEspacings, bins=nBins)

#########################
#########################

#########################
#Try an alternative method: NOT WORKING
#   Use the Gaussian Broadening Method
#########################
def gauss(x,mu,sigma):
    exponent = -1. * ((x - mu)**2) / (2.*(sigma**2))
    return denom * np.exp(exponent)

#raw_input("Paused before gaussian broadening.")
denom = (1./(sigma*np.sqrt(2.*np.pi)))

rhoGauss = np.zeros(eigvals.size)

for x1 in eigvals:
    rhoGauss = rhoGauss + gauss(eigvals,x1,sigma)

meanSpacing = np.mean(spacings)
print "Current sigma: ", sigma
print "Mean non-normalized energy spacing: ", meanSpacing
print "Mean number of states per sigma: ", sigma/meanSpacing
#########################
#########################
print "Maximum Energy: ", np.max(eigvals)
print "Minimum Energy: ", np.min(eigvals)

#Gaussian Broadening: Calculate the consistency requirement (see eqn 4.8.2 of Quantum Signatures of Chaos, 2010 by Haake):
#########################

rhoGaussDer = (np.roll(rhoGauss,-1)[1:-1] - np.roll(rhoGauss,1)[1:-1]) / ( np.abs(np.diff(eigvals)) + np.abs(np.roll(np.diff(eigvals),1)) )[1:]

#print locLevelRhoDer
checkRatioGB = np.abs(rhoGaussDer) / (rhoGauss[1:-1]**2)
#print checkRatio
print "GB: Maximum Consistency Ratio (want << 1): ", np.nanmax(checkRatioGB)

#Gaussian Broadening: Bin the data
#########################

#Normalize the spacing with the local level density
spacingsNormalizedGB = (spacings*rhoGauss[:-1])

#Bin the data
binnedDataGB, bin_edgesGB = np.histogram(spacingsNormalizedGB, bins=nBins)

print "Local level density binned."

#def subCalc():
#    rhoGaussExp = np.array([gauss(eigvals,x1,sigma) for x1 in eigvals])
#    raw_input("Paused before summation.")
#    return np.sum(rhoGaussExp,axis=0)
#
#rhoGauss = subCalc()

print "Gaussian Broadening Density calculated."
#print rhoGaussExp
#print rhoGauss
#print np.sum(rhoGaussExp[:][0])
#print rhoGauss - np.roll(rhoGauss,rhoGauss.size)

#Calculate Level density from 2nd order derivative
#O2rho = 1. / np.gradient(eigvals

print " "

#########################
#########################
#Plot
#########################
#########################
figNum = 0

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
line, = ax.plot(eigvals[nu:(-nu)], locLevelRho)#, label="Local Level Density")
#line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
#ax.legend()

#plt.ylim(0,1.5)
#plt.xlim(0,N)
plt.xlabel("Energy", fontsize=16)
plt.ylabel(r'$\rho$', fontsize=16)
plt.title("Local Level Density", fontsize=18)

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
bar = ax.bar(bin_edges[:-1], binnedData,width=bin_edges[1]-bin_edges[0])
#bar = ax.bar(bin_edges[1:-1], binnedData[1:],width=bin_edges[1]-bin_edges[0])

#plt.ylim(0,0.1)
#plt.xlim(0,20)
plt.xlabel(r'$s$', fontsize=16)
plt.ylabel(r'$P(s)$', fontsize=16)
plt.title("Local rho Energy Level Statistics", fontsize=18)

#figNum += 1
#fig = plt.figure(figNum)
#ax = plt.subplot()
#bar = ax.bar(bin_edges2[:-1], binnedData2,width=bin_edges[1]-bin_edges[0])
##bar = ax.bar(bin_edges[1:-1], binnedData[1:],width=bin_edges[1]-bin_edges[0])
#
##plt.ylim(0,0.1)
##plt.xlim(0,4)
#plt.xlabel(r'$s$', fontsize=16)
#plt.ylabel(r'$P(s)$', fontsize=16)
#plt.title("Energy Level Statistics Integrated Rho", fontsize=18)

#########################
#########################
figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
#line, = ax.plot(eigvals[nu:(-nu)], locLevelRho, label="Local Level Density")
line2, = ax.plot(eigvals, rhoGauss, color='red')#, label="Gaussian Broadening Method")
#ax.legend()

#plt.ylim(0,1.5)
#plt.xlim(0,N)
plt.xlabel("Energy", fontsize=16)
plt.ylabel(r'$\rho$', fontsize=16)
plt.title("Gaussian Broadening Level Density", fontsize=18)

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
bar = ax.bar(bin_edgesGB[:-1], binnedDataGB,width=bin_edgesGB[1]-bin_edgesGB[0],color='g',edgecolor='g')
#bar = ax.bar(bin_edges[1:-1], binnedData[1:],width=bin_edges[1]-bin_edges[0])

#plt.ylim(0,2)
#plt.xlim(0,20)
plt.xlabel(r'$s$', fontsize=16)
plt.ylabel(r'$P(s)$', fontsize=16)
plt.title("GB rho Energy Level Statistics", fontsize=18)

#########################
#########################




plt.show()



















