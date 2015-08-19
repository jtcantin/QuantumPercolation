# To run: name.py [inputDir]


import numpy as np
import sys
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import distFit
import scipy.stats as spStat

np.set_printoptions(threshold=10000,linewidth=2000,precision=4,suppress=False)

alphaCI = 0.05

inputDir = sys.argv[1]

inputFiles = [ f for f in listdir(inputDir) if (isfile(join(inputDir,f)) and f!="._.DS_Store" and f!=".DS_Store") ]



infoList = []
distArrayList = []
binEdgesArrayList = []
SValsArrayList = []
bDistArrayList = []

#Get info from each datafile
for datafileName in inputFiles:
    print "Reading: {0}".format(datafileName)
    #    datafile = open(datafileName, 'r')
    dataDict = np.load(join(inputDir,datafileName))
    
    infoList.append(dataDict['infoArray'])
    
    #    np.vstack((dist2DArray,dataDict['binnedDataGBNE']))
    distArrayList.append(dataDict['binnedDataGBNE'])
    binEdgesArrayList.append(dataDict['bin_edgesGBNE'])
#    SValsArrayList.append(dataDict['Svals'])
#    bDistArrayList.append(dataDict['brodyDistArray'])


#print infoList
distArrayArray = np.array(distArrayList)
binEdgesArrayArray = np.array(binEdgesArrayList)
#SValsArrayArray = np.array(SValsArrayList)
#bDistArrayArray = np.array(bDistArrayList)



Nlist = []
qList = []
wList = []
bValList = []
for infoArray in infoList:
    Nlist.append(int(infoArray[0]))
    qList.append(float(infoArray[1]))
    wList.append(float(infoArray[2]))
    bValList.append(float(infoArray[8]))


NArray = np.array(Nlist)
qArray = np.array(qList)
wArray = np.array(wList)
bValArray = np.array(bValList)

data2DArray = np.vstack((NArray,qArray,wArray,bValArray))

#print data2DArray

#Sort according to size and then disorder strength
#data2DArray_N_srtd_indices = np.argsort(data2DArray[2])
data2DArray_Nq_srtd_indices = np.lexsort(np.vstack((qArray,NArray)))

data2DArray_Nq_srtd = data2DArray[:,data2DArray_Nq_srtd_indices]
#print data2DArray_Nq_srtd

# Sort according to disorder strength

#data2DArray_q_srtd_indices = np.argsort(data2DArray[1])
#
#data2DArray_q_srtd = np.copy(data2DArray[:,data2DArray_q_srtd_indices])

#print data2DArray_q_srtd

#########################
#########################
#Plot
#########################
#########################
figNum = 0

#######
#b vs w
#######


figNum += 1
fig = plt.figure(figNum,facecolor="white")
ax = plt.subplot()

colourList = ["red","blue","black"]
markerList = ["D","s","."]
markerSizeList = [5,5,15]
typeIndex = 0

newNIndex = 0
NPlotList = []
bPlotList = []
errorBarList = []

while newNIndex < data2DArray_Nq_srtd[0].size:

    oldNIndex = newNIndex
    NvalueCurrent = data2DArray_Nq_srtd[0,oldNIndex]
    newNIndex = np.sum(data2DArray_Nq_srtd[0] <= data2DArray_Nq_srtd[0,newNIndex])
#    print newNIndex



    relevant_data2DArray_Nq_srtd = data2DArray_Nq_srtd[:, oldNIndex:newNIndex]
    
    NPlotList.append(relevant_data2DArray_Nq_srtd[0,0])
    
    relevant_bArray = relevant_data2DArray_Nq_srtd[3]
    
    meanB = np.mean(relevant_bArray)
    bPlotList.append(meanB)
    
    numWValues = relevant_bArray.size
    
    stdDevB = np.sqrt(np.sum((relevant_bArray - meanB)**2) / (numWValues - 1))
#   print stdDevB

    stdError = stdDevB/np.sqrt(numWValues)
    
    tDistFactor = spStat.t.ppf(1. - (alphaCI/2.), numWValues-1)
    
    errorBarHalf = tDistFactor * stdError
    
    errorBarList.append(errorBarHalf)
    
#    print relevant_data2DArray_Nq_srtd

NPlotArray = np.array(NPlotList)
bPlotArray = np.array(bPlotList)
errorBarArray = np.array(errorBarList)





#print data2DArray_Nq_srtd[0] <= data2DArray_Nq_srtd[0,0]
#print np.sum(data2DArray_Nq_srtd[0] <= data2DArray_Nq_srtd[0,0]) #This gives next index after first value
#print data2DArray_Nq_srtd[:,data2DArray_Nq_srtd[0] <= data2DArray_Nq_srtd[0,0]]







#yerr=errorBarArray
line = ax.errorbar(NPlotArray, bPlotArray, yerr=errorBarArray, marker=markerList[typeIndex], markersize=markerSizeList[typeIndex], label="N={0}".format(NvalueCurrent), clip_on=False, linewidth=1.0, color=colourList[typeIndex])#linewidth=3.0, ls=''
typeIndex += 1
#    line[-1][0].set_linestyle('--')
#    line, = ax.errorbar(qPlotList, bPlotList, yerr=errorBarList, 'k-',marker=".", markersize=20,linewidth=3.0, label="Diffusive Regime", clip_on=False)
#line, = ax.plot(data2DArray_q_srtd[1], data2DArray_q_srtd[3],'k-',marker=".", markersize=20,linewidth=3.0, label="Diffusive Regime", clip_on=False)
#line2, = ax.plot(second_q_NArray, second_q_bValArray,'k-',marker=".",markersize=20,linewidth=2.0, label="Localized Regime", clip_on=False)
#line, = ax.plot(NArray_srtd, bValArray_srtd,'k-',marker=".",markersize=20,linewidth=2.0)#, label="Local Level Density")
#line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
#ax.legend()

# Remove plot frame
#ax.set_frame_on(False)



#plt.ylim(0,1.)
plt.xlim(25,85)
plt.xlabel("System Side Length, N", fontsize=16)
plt.ylabel("Brody Parameter, b", fontsize=16)
#plt.title("Local Level Density", fontsize=18)

print "N:"
print NPlotArray
print "Mean b:"
print bPlotArray

plt.show()




