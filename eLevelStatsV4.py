# To run: python eLevelStatsV2.py [eigvalFilename] [nu] [nBins] [sigma] [sigma2] 

import sys
import numpy as np
from scipy import integrate, optimize
import matplotlib.pyplot as plt
import time
import distFit

np.set_printoptions(threshold=10000,linewidth=2000,precision=20,suppress=False)

ti = time.time()

eigvalFilename = sys.argv[1]
nu = int(sys.argv[2])
nBins = int(sys.argv[3])
sigma = float(sys.argv[4])
sigma2 = float(sys.argv[5])
#sigma3 = float(sys.argv[6])
#bVal = float(sys.argv[6])
#NVal = int(sys.argv[6])
#qVal = float(sys.argv[7])
#wVal = float(sys.argv[8])
#rangeString = sys.argv[9]

plotCheck = True

###########
#Get parameters from directory name
###########
locNameSplit = eigvalFilename.split('/')
dirName = locNameSplit[-2]
NVal = int(dirName[6:9])
qVal = float(dirName[11:13])/100.
wVal = float(dirName[15:18])
rangeString = str(dirName[20:22])
dateString = str(dirName[23:])

print "N: ", NVal
print "q: ", qVal
print "w: ", wVal
print "Range: ", rangeString
print "Date: ", dateString
junk = raw_input("Check these. Hit Enter to continue.")



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
binnedData, bin_edges = np.histogram(spacingsNormalized, bins=nBins,density=True)

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
#   Calculate staircase function
#########################
def staircase(x,eigvals):
    total = x.size
    NofEreal = np.zeros(total)

    for i in range(0,total):
        NofEreal[i] = np.sum(eigvals<(x[i]))

    return NofEreal

EValues = np.arange(np.min(eigvals), np.max(eigvals),0.0001)
staircaseArray = staircase(EValues, eigvals)

#########################
#Try an alternative method: NOT WORKING
#   Use the Gaussian Broadening Method
#########################
def gauss(x,mu,sigma,denom):
    exponent = -1. * ((x - mu)**2) / (2.*(sigma**2))
    return denom * np.exp(exponent)

#raw_input("Paused before gaussian broadening.")
#denom = (1./(sigma*np.sqrt(2.*np.pi)))

def gaussSum2(x,eigvals,sigma):
    rhoGauss = np.zeros(len(x))
    denom = (1./(sigma*np.sqrt(2.*np.pi)))
    
    for x1 in eigvals: #Sum over the different gaussians, each with a mean equal to the eigenvalue
        rhoGauss = rhoGauss + gauss(x,x1,sigma,denom)

    return rhoGauss

def gaussSum(x,eigvals,sigma):
    rhoGauss = 0.
    denom = (1./(sigma*np.sqrt(2.*np.pi)))
    
    for x1 in eigvals: #Sum over the different gaussians, each with a mean equal to the eigenvalue
        rhoGauss = rhoGauss + gauss(x,x1,sigma,denom)
    
    return rhoGauss
#Evals = np.arange(np.min(eigvals)-1, np.max(eigvals)+1,0.0001)
#print Evals.size
Evals = eigvals

rhoGauss = gaussSum2(Evals,eigvals,sigma)

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

rhoGaussDer = (np.roll(rhoGauss,-1)[1:-1] - np.roll(rhoGauss,1)[1:-1]) / ( np.abs(np.diff(Evals)) + np.abs(np.roll(np.diff(Evals),1)) )[1:]

#print locLevelRhoDer
checkRatioGB = np.abs(rhoGaussDer) / (rhoGauss[1:-1]**2)
#print checkRatio
print "GB: Maximum Consistency Ratio (want << 1): ", np.nanmax(checkRatioGB)

#Gaussian Broadening: Bin the data
#########################

#Normalize the spacing with the local level density
spacingsNormalizedGB = (spacings*rhoGauss[:-1])



#Bin the data
binnedDataGB, bin_edgesGB = np.histogram(spacingsNormalizedGB, bins=nBins,density=True)

print "Gaussian Broadened density binned."

#def subCalc():
#    rhoGaussExp = np.array([gauss(eigvals,x1,sigma) for x1 in eigvals])
#    raw_input("Paused before summation.")
#    return np.sum(rhoGaussExp,axis=0)
#
#rhoGauss = subCalc()

print "Gaussian Broadening Density calculated."

#GB: Integrate to get mean N(E)
###############
print "Integrating Gaussian Broadened rho."

NoFE_array = integrate.cumtrapz(rhoGauss,Evals,initial=0.5)

print "Integration Complete."

def MSD_NofE(sigma, Evals, eigvals, testArray, arrayReturn=False):
    rhoGauss = gaussSum2(Evals,eigvals,sigma)
    NoFE_array = integrate.cumtrapz(rhoGauss,Evals,initial=0.5)

    Dev = NoFE_array - testArray
    SqrDev = Dev**2
    meanSqrDev = np.mean(SqrDev)

    if arrayReturn:
        return meanSqrDev, NoFE_array, rhoGauss

    else:
        return meanSqrDev

#GB: Use the secant method to obtain the least mean squares from the staircase function
###############
print "Time so far: ", time.time() - ti, " s"
ti = time.time()

Evals = eigvals
staircaseArrayTest = staircase(Evals, eigvals)

#Minimize MSD
print "Beginning minimization..."
ti = time.time()

MSDoptions = {"maxiter": 3, "disp": True}
MinTol = 1E-10 #Relative tolerance
print "Chosen Relative Tolerance: ", MinTol

MSDresult = optimize.minimize_scalar(MSD_NofE, bracket=[sigma, sigma2], args=(Evals, eigvals, staircaseArrayTest), tol=MinTol)#, options=MSDoptions)

print "Optimization Complete. Time: ", time.time() - ti, " s"
ti = time.time()
print "MSDresult:"
print MSDresult
print " "
#Original sigma
Dev = NoFE_array - staircaseArrayTest
SqrDev = Dev**2
meanSqrDev = np.mean(SqrDev)

sigmaFinal = MSDresult.x

print "Original Sigma: ", sigma
print "Original MSD: ", meanSqrDev
print "Determined Sigma: ", MSDresult.x
#print "Optimizer Exited Successfully: ", MSDresult.success
#print "Obtained MSD: ", MSDresult.fun
#print "Number of iterations: ", MSDresult.nit

#
#
##Second Sigma
#meanSqrDev2, NoFE_array2 = MSD_NofE(sigma2, Evals, eigvals, staircaseArrayTest)
##rhoGauss2 = gaussSum2(Evals,eigvals,sigma2)
##NoFE_array2 = integrate.cumtrapz(rhoGauss2,Evals,initial=0.5)
##
##Dev = NoFE_array2 - staircaseArrayTest
##SqrDev = Dev**2
##meanSqrDev2 = np.mean(SqrDev)
#
#
##Third Sigma
#meanSqrDev3, NoFE_array2 = MSD_NofE(sigma3, Evals, eigvals, staircaseArrayTest)
##rhoGauss3 = gaussSum2(Evals,eigvals,sigma3)
##NoFE_array3 = integrate.cumtrapz(rhoGauss3,Evals,initial=0.5)
##
##Dev = NoFE_array3 - staircaseArrayTest
##SqrDev = Dev**2
##meanSqrDev3 = np.mean(SqrDev)
#
#
##Calculate Derivatives
#Der1 = (meanSqrDev2 - meanSqrDev) / (sigma2 - sigma)
#Der2 = (meanSqrDev3 - meanSqrDev2) / (sigma3 - sigma2)
#
#print "Time for prep: ", time.time() - ti, "s"
#ti = time.time()
#print "Sigma 1: ", sigma, " MSD: ", meanSqrDev
#print "Sigma 2: ", sigma2, " MSD: ", meanSqrDev2, "Derivative: ", Der1
#print "Sigma 3: ", sigma3, " MSD: ", meanSqrDev3, "Derivative: ", Der2
#
#figNum = 0
#
##Use the secant method until relative error achieved
#relError = 0.01
##while np.abs((meanSqrDev2 - meanSqrDev)/meanSqrDev) > relError:
#while np.abs((sigma3 - sigma2)/sigma2) > relError:
#    
#    #Get new sigma
#    sigmaNew = sigma3 - Der2 * (sigma3 - sigma2) / (Der2 - Der1)
#
#    Der1 = Der2
#
#    meanSqrDev2 = meanSqrDev3
#
#    #Calculate with new sigma
#    meanSqrDev3, NoFE_array2 = MSD_NofE(sigmaNew, Evals, eigvals, staircaseArrayTest)
##    rhoGauss2 = gaussSum2(Evals,eigvals,sigmaNew)
##    NoFE_array2 = integrate.cumtrapz(rhoGauss2,Evals,initial=0.5)
##
##    #Get new mean squared deviation
##    Dev = NoFE_array2 - staircaseArrayTest
##    SqrDev = Dev**2
##    meanSqrDev3 = np.mean(SqrDev)
#
#    Der2 = (meanSqrDev3 - meanSqrDev2) / (sigmaNew - sigma3)
#
#    #Update sigmas
#    sigma2 = sigma3
#    sigma3 = sigmaNew
#
#    print "New Sigma: ", sigmaNew, " New MSD: ", meanSqrDev3, "Derivative: ", Der2, " Iteration Time: ", time.time() - ti, "s"
#    
#    if plotCheck:
#        figNum += 1
#        fig = plt.figure(figNum)
#        ax = plt.subplot()
#        line, = ax.plot(Evals, NoFE_array2, color='blue', label="meanN(E)")
#        line2, = ax.plot(EValues, staircaseArray, color='red', label="N(E)")#, label="Gaussian Broadening Method")
#        ax.legend()
#
#        #plt.ylim(0,1.5)
#        #plt.xlim(0,N)
#        plt.xlabel("Energy", fontsize=16)
#        plt.ylabel(r'$N(E)$', fontsize=16)
#        plt.title("Integrated Gaussian Broadening Level Density", fontsize=18)
#        
#        plt.show()
##        raw_input("Check out plot.")
#
#    
#
#    ti = time.time()
#
#NoFE_array = NoFE_array2

meanSqrDevFinal, NoFE_array, rhoGauss2 = MSD_NofE(MSDresult.x, Evals, eigvals, staircaseArrayTest, arrayReturn=True)

sigmaRMSD = np.sqrt(meanSqrDevFinal)

#Gaussian Broadening: Bin the unfolded data
#########################

#Normalize the spacing with the local level density
#spacingsNormalizedGB = (spacings*rhoGauss[:-1])

spacingsNormalizedGB = np.diff(NoFE_array)

#Check for any negative values
print "The following negative values encountered in NoFE_GB:"
print spacingsNormalizedGB[spacingsNormalizedGB < 0]
print "Locations: "
print np.where(spacingsNormalizedGB < 0)


if (spacingsNormalizedGB[spacingsNormalizedGB < 0]).size > 0:
    print NoFE_array[0]
    print NoFE_array[1]
    raw_input("ERROR: Negative spacing")

#Bin the data
binnedDataGBNE, bin_edgesGBNE = np.histogram(spacingsNormalizedGB, bins=nBins,density=True)

print "Gaussian Broadened staricase binned."

#total = eigvals.size
#lowerBound = rhoGauss
#NoFE_array = np.zeros(total)
#
#NofE, abserr = integrate.quad(gaussSum, np.min(eigvals)-0.01, eigvals[0], args=(eigvals,sigma))
#NoFE_array[0] = NofE
#print "%d of %d integrations completed." % (1, total)
#
#for i in range(1,total):
#    NofE, abserr = integrate.quad(gaussSum, eigvals[i-1], eigvals[i], args=(eigvals,sigma))
#    NoFE_array[i] = NoFE_array[i-1] + NofE
#    print "%d of %d integrations completed." % (i+1, total)


#NofE, abserr = integrate.quad(gaussSum, np.min(eigvals)-0.01, 0., args=(eigvals,sigma))

#print NofE
#print abserr

#print rhoGaussExp
#print rhoGauss
#print np.sum(rhoGaussExp[:][0])
#print rhoGauss - np.roll(rhoGauss,rhoGauss.size)

#Calculate Level density from 2nd order derivative
#O2rho = 1. / np.gradient(eigvals

#Gaussian Broadening: Fit to Brody Distribution
#########################

#Set up the points to sit at the centre of each bin
Svals = bin_edgesGBNE[:-1] + (bin_edgesGBNE[1]-bin_edgesGBNE[0])/2.

#brodyDistArray = distFit.brodyDist(Svals, bVal)

def MSD_BrodyDist(bVal, Svals, dataToFitArray, arrayReturn=False):
    brodyDistArray = distFit.brodyDist(Svals, bVal)

    devArray = dataToFitArray - brodyDistArray
    sqrDevArray = np.power(devArray, 2)
    meanSqrDev = np.mean(sqrDevArray)

    if arrayReturn:
        return meanSqrDev, brodyDistArray
    
    else:
        return meanSqrDev

#Fit Brody Dist by minimizing MSD
print "Beginning Brody Distribution Fit..."
ti = time.time()

MSDoptions = {"maxiter": 3, "disp": True}
MinTol = 1E-10 #Relative tolerance
print "Chosen Relative Tolerance: ", MinTol

fitResult = optimize.minimize_scalar(MSD_BrodyDist, bracket=[0., 1.], args=(Svals, binnedDataGBNE), tol=MinTol)#, options=MSDoptions)

print "Fit Complete. Time: ", time.time() - ti, " s"
ti = time.time()

bVal = fitResult.x

print "Fit Result:"
print fitResult
print "Determined b value: ", bVal

MSD_fit, brodyDistArray = MSD_BrodyDist(bVal, Svals, binnedDataGBNE, arrayReturn=True)
print "MSD of the Fit: ", MSD_fit

RMSD_fit = np.sqrt(MSD_fit)
bValRMSD = RMSD_fit

print " "

#########################
#########################
#Save Data
#########################
#########################
infoArray = np.array([NVal,qVal,wVal,rangeString,nBins,nu,sigmaFinal,sigmaRMSD,bVal,bValRMSD])
infoArrayName = "[NVal,qVal,wVal,rangeString,nBins,nu,sigmaFinal,sigmaRMSD,bVal,bValRMSD]"

distArrayFilename = eigvalFilename[:-4] + "_distArrayData.npz"

np.savez(distArrayFilename, infoArray="[NVal,qVal,wVal,rangeString,nBins,nu,sigmaFinal,bVal]", binnedDataGBNE=binnedDataGBNE, bin_edgesGBNE=bin_edgesGBNE, brodyDistArray=brodyDistArray, Svals=Svals)

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
#line, = ax.plot(Evals[1:], NoFE_array, color='blue', label="N(E)")
line2, = ax.plot(Evals, rhoGauss, color='red', label="Original Sigma")
line2, = ax.plot(Evals, rhoGauss2, color='blue', label="Final Sigma")
ax.legend()

#plt.ylim(0,1.5)
#plt.xlim(0,N)
plt.xlabel("Energy", fontsize=16)
plt.ylabel(r'$\rho$', fontsize=16)
plt.title("Gaussian Broadening Level Density", fontsize=18)

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
line, = ax.plot(Evals, NoFE_array, color='blue', label="meanN(E)")
line2, = ax.plot(EValues, staircaseArray, color='red', label="N(E)")#, label="Gaussian Broadening Method")
ax.legend()

#plt.ylim(0,1.5)
#plt.xlim(0,N)
plt.xlabel("Energy", fontsize=16)
plt.ylabel(r'$N(E)$', fontsize=16)
plt.title("Integrated Gaussian Broadening Level Density", fontsize=18)

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
#bar = ax.bar(bin_edgesGB[:-1], binnedDataGB,width=bin_edgesGB[1]-bin_edgesGB[0],color='g',edgecolor='g',label="rho unfolded")
#bar2 = ax.bar(bin_edgesGBNE[:-1], binnedDataGBNE,width=bin_edgesGBNE[1]-bin_edgesGBNE[0],color='k',edgecolor='k',label="N(E_k)")
line2, = ax.plot(bin_edgesGBNE[:-1], binnedDataGBNE, 'ok', markersize=4, label="N(E_k)", markerfacecolor='w')


brodyLabel = "Brody Dist, b = %f" % bVal
line, = ax.plot(Svals, brodyDistArray, '.r', label=brodyLabel, markersize=2)
#bar = ax.bar(bin_edges[1:-1], binnedData[1:],width=bin_edges[1]-bin_edges[0])
ax.legend()

titleString = "GB Energy Level Statistics for N=%02d,q=%.2f,w=%.2f %s; RMSD=%.2f" % (NVal, qVal, wVal, rangeString, RMSD_fit)

plt.ylim(0,1.2)
#plt.xlim(0,20)
plt.xlabel(r'$s$', fontsize=16)
plt.ylabel(r'$P(s)$', fontsize=16)
plt.title(titleString, fontsize=18)



#########################
#########################




plt.show()



















