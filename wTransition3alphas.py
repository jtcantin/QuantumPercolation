# To run: name.py [inputDir]


import numpy as np
import sys
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import distFit
import scipy.stats as spStat
import scipy.interpolate

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
rangeList = []
for infoArray in infoList:
    Nlist.append(int(infoArray[0]))
    qList.append(float(infoArray[1]))
    wList.append(float(infoArray[2]))
    bValList.append(float(infoArray[8]))
    
    rangeString = infoArray[3].strip("-")
    if (rangeString == "rT") or (rangeString == "TB"):
        rangeList.append(np.inf)
    else:
        rangeList.append(float(infoArray[3].strip("-")))


NArray = np.array(Nlist)
qArray = np.array(qList)
wArray = np.array(wList)
bValArray = np.array(bValList)
rangeArray = np.array(rangeList)

data2DArray = np.vstack((NArray,qArray,wArray,bValArray,rangeArray))

#print data2DArray

#Sort according to range and then disorder strength
data2DArray_rw_srtd_indices = np.lexsort(np.vstack((wArray,rangeArray)))

data2DArray_rw_srtd = data2DArray[:,data2DArray_rw_srtd_indices]
#print data2DArray_rw_srtd



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
markerList2 = ["D","s",".","^"]
markerList = ["None","None","None","None"]
LineStyleList = [":","-.","--","-"]
markerSizeList = [5,5,15,10]
#markerSizeList = [2,2,5,5]
typeIndex = 0

newrIndex = 0
relValFirstIndex = 4

while newrIndex < data2DArray_rw_srtd[relValFirstIndex].size:

    oldrIndex = newrIndex
    rvalueCurrent = data2DArray_rw_srtd[relValFirstIndex,oldrIndex]
    newrIndex = np.sum(data2DArray_rw_srtd[relValFirstIndex] <= data2DArray_rw_srtd[relValFirstIndex,newrIndex])
#    print newrIndex

    relevant_data2DArray_rw_srtd = data2DArray_rw_srtd[:, oldrIndex:newrIndex]
    
#    print relevant_data2DArray_rw_srtd

    newWIndex = 0
    wPlotList = []
    bPlotList = []
    errorBarList = []
    stdDevBList = []

    while newWIndex < relevant_data2DArray_rw_srtd[2].size:

        oldWIndex = newWIndex
        wPlotList.append(relevant_data2DArray_rw_srtd[2,oldWIndex])
        
        newWIndex += np.sum(np.abs(relevant_data2DArray_rw_srtd[2] - relevant_data2DArray_rw_srtd[2,newWIndex]) < 1E-10)
#        print newWIndex

        relevant_bArray = relevant_data2DArray_rw_srtd[3,oldWIndex:newWIndex]
#        print relevant_bArray

        meanB = np.mean(relevant_bArray)
        bPlotList.append(meanB)
#        print meanB


        numWValues = relevant_bArray.size

        if numWValues < 2: #TEMPORARY FIX!
            numWValues +=1

        stdDevB = np.sqrt(np.sum((relevant_bArray - meanB)**2) / (numWValues - 1))
        stdDevBList.append(stdDevB)
#        print stdDevB

        stdError = stdDevB/np.sqrt(numWValues)
        
        tDistFactor = spStat.t.ppf(1. - (alphaCI/2.), numWValues-1)

        errorBarHalf = tDistFactor * stdError

        errorBarList.append(errorBarHalf)

#    print wPlotList
#    print bPlotList
#    print errorBarList
    wPlotArray = np.array(wPlotList)
    bPlotArray = np.array(bPlotList)
    errorBarArray = np.array(errorBarList)
    stdDevBArray = np.array(stdDevBList)

#print data2DArray_rw_srtd[0] <= data2DArray_rw_srtd[0,0]
#print np.sum(data2DArray_rw_srtd[0] <= data2DArray_rw_srtd[0,0]) #This gives next index after first value
#print data2DArray_rw_srtd[:,data2DArray_rw_srtd[0] <= data2DArray_rw_srtd[0,0]]







#yerr=errorBarArray
#    line = ax.errorbar(wPlotArray, bPlotArray, yerr=errorBarArray, marker=markerList[typeIndex], markersize=markerSizeList[typeIndex], label=r'$\alpha={0}$'.format(rvalueCurrent), clip_on=False, linewidth=1.0, color=colourList[typeIndex])#linewidth=3.0, ls=''

    if rvalueCurrent == np.inf:
        labelText = r'$NN$'
    else:
        labelText = r'$\alpha={0}$'.format(rvalueCurrent)

    #Spline Interpolation
    splineXarray = np.linspace(wPlotArray.min(), wPlotArray.max(), num=1000)
    akimaSplineInterpFunction = scipy.interpolate.Akima1DInterpolator(wPlotArray, bPlotArray)
#    splineInterpFunction = scipy.interpolate.interp1d(wPlotArray, bPlotArray, kind='cubic')

#    splineYarray = splineInterpFunction(splineXarray)
    splineYarray = akimaSplineInterpFunction(splineXarray)


#    tck = scipy.interpolate.splrep(wPlotArray, bPlotArray, w=(1./stdDevBArray), k=3)
#    splineYarray = scipy.interpolate.splev(splineXarray, tck, der=0)
#    print wPlotArray.size-np.sqrt(2.*wPlotArray.size)
#    print wPlotArray.size+np.sqrt(2.*wPlotArray.size)

    print "Alpha = ", rvalueCurrent
    print wPlotArray
    print bPlotArray
#    print splineXarray.T
    dataToSave = np.c_[splineXarray, splineYarray]


    np.savetxt("splineData_alpha{0}.csv".format(rvalueCurrent), dataToSave, delimiter=",", header="Akima Spline interpolation data for alpha = {0}. The first column is w, the second is b.".format(rvalueCurrent))

#    line = ax.plot(splineXarray, splineYarray, marker=markerList[typeIndex], markersize=markerSizeList[typeIndex], label=labelText, clip_on=False, linewidth=3.0, color=colourList[typeIndex], linestyle=LineStyleList[typeIndex])

    line = ax.plot(splineXarray, splineYarray, marker="None", markersize=markerSizeList[typeIndex], label=labelText, clip_on=False, linewidth=3.0, color=colourList[typeIndex], linestyle=LineStyleList[typeIndex])

    #Original Plot
#    line = ax.plot(wPlotArray, bPlotArray, marker=markerList[typeIndex], markersize=markerSizeList[typeIndex], label=labelText, clip_on=False, linewidth=3.0, color=colourList[typeIndex], linestyle=LineStyleList[typeIndex])#linewidth=3.0, ls=''
#    line2 = ax.plot(wPlotArray, bPlotArray, marker=markerList2[typeIndex], markersize=markerSizeList[typeIndex], label=labelText, clip_on=False, linewidth=3.0, color=colourList[typeIndex], linestyle='None')#linewidth=3.0, ls=''

    typeIndex += 1

#    line[-1][0].set_linestyle('--')
#    line, = ax.errorbar(wPlotList, bPlotList, yerr=errorBarList, 'k-',marker=".", markersize=20,linewidth=3.0, label="Diffusive Regime", clip_on=False)
#line, = ax.plot(data2DArray_w_srtd[2], data2DArray_w_srtd[3],'k-',marker=".", markersize=20,linewidth=3.0, label="Diffusive Regime", clip_on=False)
#line2, = ax.plot(second_q_NArray, second_q_bValArray,'k-',marker=".",markersize=20,linewidth=2.0, label="Localized Regime", clip_on=False)
#line, = ax.plot(NArray_srtd, bValArray_srtd,'k-',marker=".",markersize=20,linewidth=2.0)#, label="Local Level Density")
#line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
ax.legend()
legend = ax.legend()
legend.get_frame().set_linewidth(2.0)

# Remove plot frame
#ax.set_frame_on(False)



plt.ylim(0,1.)
#plt.xlim(0,N)
#plt.xlabel("Disorder Strength, w", fontsize=16)
#plt.ylabel("Brody Parameter, b", fontsize=16)
plt.xlabel(r'$w$', fontsize=24)
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

figureFilename = "wTransitionAlphasN25_q25.eps"
#fig.savefig(figureFilename, format='eps', dpi=1200)


plt.show()




