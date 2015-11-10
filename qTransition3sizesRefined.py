# To run: name.py [inputDir]


import numpy as np
import sys
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import distFit
import scipy.stats as spStat
import scipy.interpolate
import matplotlib.ticker

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

#colourList = ["black","blue","green","red"]
colourList = ["black","black","black","black"]
markerList = ["D","s",".","^"]
#markerList = ["None","None","None","None"]
LineStyleList = [":","--","-","-."]
markerSizeList = [5,5,15,10]
#markerSizeList = [2,2,5,5]
markerSizeList = [5,5,15]
typeIndex = 0

newNIndex = 0

while newNIndex < data2DArray_Nq_srtd[0].size:

    oldNIndex = newNIndex
    NvalueCurrent = data2DArray_Nq_srtd[0,oldNIndex]
    newNIndex = np.sum(data2DArray_Nq_srtd[0] <= data2DArray_Nq_srtd[0,newNIndex])
#    print newNIndex

    relevant_data2DArray_Nq_srtd = data2DArray_Nq_srtd[:, oldNIndex:newNIndex]
    
#    print relevant_data2DArray_Nq_srtd

    newWIndex = 0
    qPlotList = []
    bPlotList = []
    errorBarList = []

    while newWIndex < relevant_data2DArray_Nq_srtd[1].size:

        oldWIndex = newWIndex
        qPlotList.append(relevant_data2DArray_Nq_srtd[1,oldWIndex])
        
        newWIndex += np.sum(np.abs(relevant_data2DArray_Nq_srtd[1] - relevant_data2DArray_Nq_srtd[1,newWIndex]) < 1E-10)
#        print newWIndex

        relevant_bArray = relevant_data2DArray_Nq_srtd[3,oldWIndex:newWIndex]
#        print relevant_bArray

        meanB = np.mean(relevant_bArray)
        bPlotList.append(meanB)
#        print meanB


        numWValues = relevant_bArray.size

#        if numWValues < 2: #TEMPORARY FIX!
#            numWValues +=1

        stdDevB = np.sqrt(np.sum((relevant_bArray - meanB)**2) / (numWValues - 1))
#        print stdDevB

        stdError = stdDevB/np.sqrt(numWValues)
        
        tDistFactor = spStat.t.ppf(1. - (alphaCI/2.), numWValues-1)

        errorBarHalf = tDistFactor * stdError

        errorBarList.append(errorBarHalf)

#    print qPlotList
#    print bPlotList
#    print errorBarList
    qPlotArray = np.array(qPlotList)
    pPlotArray = 1.-qPlotArray[::-1].copy()
    bPlotArray = np.array(bPlotList)
    bPlotArray_rev = bPlotArray[::-1].copy()
    errorBarArray = np.array(errorBarList)

#Akima Spline Interpolation
    splineXarray = np.linspace(pPlotArray.min(), pPlotArray.max(), num=1000)
    akimaSplineInterpFunction = scipy.interpolate.Akima1DInterpolator(pPlotArray, bPlotArray_rev)
    
    splineYarray = akimaSplineInterpFunction(splineXarray)

#print data2DArray_Nq_srtd[0] <= data2DArray_Nq_srtd[0,0]
#print np.sum(data2DArray_Nq_srtd[0] <= data2DArray_Nq_srtd[0,0]) #This gives next index after first value
#print data2DArray_Nq_srtd[:,data2DArray_Nq_srtd[0] <= data2DArray_Nq_srtd[0,0]]







#yerr=errorBarArray
#    line = ax.errorbar(qPlotArray, bPlotArray, yerr=errorBarArray, marker=markerList[typeIndex], markersize=markerSizeList[typeIndex], label="N={0}".format(NvalueCurrent), clip_on=False, linewidth=1.0, color=colourList[typeIndex])#linewidth=3.0, ls=''

    line = ax.plot(splineXarray, splineYarray, marker="None", markersize=markerSizeList[typeIndex], label=r'$N={0:d}$'.format(int(NvalueCurrent)), clip_on=False, linewidth=3.0, color=colourList[typeIndex], linestyle=LineStyleList[typeIndex])#linewidth=3.0, ls=''

    typeIndex += 1
#    line[-1][0].set_linestyle('--')
#    line, = ax.errorbar(qPlotList, bPlotList, yerr=errorBarList, 'k-',marker=".", markersize=20,linewidth=3.0, label="Diffusive Regime", clip_on=False)
#line, = ax.plot(data2DArray_q_srtd[1], data2DArray_q_srtd[3],'k-',marker=".", markersize=20,linewidth=3.0, label="Diffusive Regime", clip_on=False)
#line2, = ax.plot(second_q_NArray, second_q_bValArray,'k-',marker=".",markersize=20,linewidth=2.0, label="Localized Regime", clip_on=False)
#line, = ax.plot(NArray_srtd, bValArray_srtd,'k-',marker=".",markersize=20,linewidth=2.0)#, label="Local Level Density")
#line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
ax.legend()
legend = ax.legend(loc='upper left')
legend.get_frame().set_linewidth(2.0)

# Remove plot frame
#ax.set_frame_on(False)

majorLocator = matplotlib.ticker.MultipleLocator(0.2)
ax.xaxis.set_major_locator(majorLocator)

plt.ylim(0,1.)
plt.xlim(0,1.)
#plt.xlabel("Dilution, q", fontsize=16)
#plt.ylabel("Brody Parameter, b", fontsize=16)
plt.xlabel(r'$p$', fontsize=24)
plt.ylabel(r'$b$', fontsize=24)
#plt.title("Local Level Density", fontsize=18)
[i.set_linewidth(3.0) for i in ax.spines.itervalues()]
for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(6.)
    tick.label1 = tick._get_text1()
for tick in ax.get_yaxis().get_major_ticks():
    tick.set_pad(8.)
    tick.label1 = tick._get_text1()
for label in ax.yaxis.get_ticklabels():
    label.set_verticalalignment('center')
ax.tick_params(direction="inout", length=10, width=2, colors='k', top='off', right='off', labelsize=20)
plt.tight_layout()

figureFilename = "/Users/jtcantin/Documents/qTransitionNN_w004.eps"
fig.savefig(figureFilename, format='eps', dpi=1200)

plt.show()




