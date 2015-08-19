# To run: name.py [inputDir]


import numpy as np
import sys
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import distFit

np.set_printoptions(threshold=10000,linewidth=2000,precision=4,suppress=False)

inputDir = sys.argv[1]

inputFiles = [ f for f in listdir(inputDir) if (isfile(join(inputDir,f)) and f!="._.DS_Store" and f!=".DS_Store") ]

infoList = []
distArrayList = []
binEdgesArrayList = []
SValsArrayList = []
bDistArrayList = []

#Get info from each datafile
for datafileName in inputFiles:
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

# Sort according to disorder strength

data2DArray_w_srtd_indices = np.argsort(data2DArray[2])

data2DArray_w_srtd = np.copy(data2DArray[:,data2DArray_w_srtd_indices])

#print data2DArray_w_srtd

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
line, = ax.plot(data2DArray_w_srtd[2], data2DArray_w_srtd[3],'k-',marker=".", markersize=20,linewidth=3.0, label="Diffusive Regime", clip_on=False)
#line2, = ax.plot(second_q_NArray, second_q_bValArray,'k-',marker=".",markersize=20,linewidth=2.0, label="Localized Regime", clip_on=False)
#line, = ax.plot(NArray_srtd, bValArray_srtd,'k-',marker=".",markersize=20,linewidth=2.0)#, label="Local Level Density")
#line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
#ax.legend()

# Remove plot frame
#ax.set_frame_on(False)



plt.ylim(0,1.)
#plt.xlim(0,N)
plt.xlabel("Disorder Strength, w", fontsize=16)
plt.ylabel("Brody Parameter, b", fontsize=16)
#plt.title("Local Level Density", fontsize=18)

plt.show()




