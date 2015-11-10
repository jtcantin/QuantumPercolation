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
data2DArray_Nqw_srtd_indices = np.lexsort(np.vstack((wArray,qArray,NArray)))

data2DArray_Nqw_srtd = data2DArray[:,data2DArray_Nqw_srtd_indices]
#print data2DArray_Nqw_srtd

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

#cm = plt.cm.get_cmap('RdYlBu')
#cm = plt.cm.get_cmap('YlGnBu_r')
#cm = plt.cm.get_cmap('YlOrBr_r')
cm = plt.cm.get_cmap('hot')
#cm = plt.cm.get_cmap('gray')
#cm = plt.cm.get_cmap('')

colourList = ["red","blue","black"]
markerList = ["D","s","."]
markerSizeList = [5,5,15]
typeIndex = 0

newNIndex = 0

plotHeight = 21
plotWidth = 1.

while newNIndex < data2DArray_Nqw_srtd[0].size:

    oldNIndex = newNIndex
    NvalueCurrent = data2DArray_Nqw_srtd[0,oldNIndex]
    newNIndex = np.sum(data2DArray_Nqw_srtd[0] <= data2DArray_Nqw_srtd[0,newNIndex])
#    print newNIndex

    relevant_data2DArray_Nqw_srtd = data2DArray_Nqw_srtd[:, oldNIndex:newNIndex]
    
#    print relevant_data2DArray_Nqw_srtd


    newQIndex = 0
    wPlotList = []
    qPlotList = []
    bPlotList = []
    errorBarList = []
    
    while newQIndex < relevant_data2DArray_Nqw_srtd[1].size:
        oldQIndex = newQIndex
#        qPlotList.append(relevant_data2DArray_Nqw_srtd[1,oldQIndex])

        newQIndex += np.sum(np.abs(relevant_data2DArray_Nqw_srtd[1] - relevant_data2DArray_Nqw_srtd[1,newQIndex]) < 1E-10)
    
        w_relevant_data2DArray_Nqw_srtd = relevant_data2DArray_Nqw_srtd[:, oldQIndex:newQIndex]
    
        newWIndex = 0

        while newWIndex < w_relevant_data2DArray_Nqw_srtd[2].size:

            oldWIndex = newWIndex
            wPlotList.append(w_relevant_data2DArray_Nqw_srtd[2,oldWIndex])
            qPlotList.append(relevant_data2DArray_Nqw_srtd[1,oldQIndex])

            newWIndex += np.sum(np.abs(w_relevant_data2DArray_Nqw_srtd[2] - w_relevant_data2DArray_Nqw_srtd[2,newWIndex]) < 1E-10)

            relevant_bArray = w_relevant_data2DArray_Nqw_srtd[3,oldWIndex:newWIndex]

            meanB = np.mean(relevant_bArray)
            bPlotList.append(meanB)

            numDataValues = relevant_bArray.size
            stdDevB = np.sqrt(np.sum((relevant_bArray - meanB)**2) / (numDataValues - 1))
            stdError = stdDevB/np.sqrt(numDataValues)

            tDistFactor = spStat.t.ppf(1. - (alphaCI/2.), numDataValues-1)
                    
            errorBarHalf = tDistFactor * stdError
                        
            errorBarList.append(errorBarHalf)

#    print wPlotList
#    print qPlotList
#    print bPlotList
#    print errorBarList
    wPlotArray = np.array(wPlotList)
    qPlotArray = np.array(qPlotList)
    bPlotArray = np.array(bPlotList)
    errorBarArray = np.array(errorBarList)

#print data2DArray_Nqw_srtd[0] <= data2DArray_Nqw_srtd[0,0]
#print np.sum(data2DArray_Nqw_srtd[0] <= data2DArray_Nqw_srtd[0,0]) #This gives next index after first value
#print data2DArray_Nqw_srtd[:,data2DArray_Nqw_srtd[0] <= data2DArray_Nqw_srtd[0,0]]






    print qPlotArray
#    print qPlotArray[::-1]
    print wPlotArray
    print bPlotArray
#    bPlot2DArrayOld = np.reshape(np.append(bPlotArray,np.nan),(5,-1)) #r5
    bPlot2DArrayOld = np.reshape(bPlotArray,(4,-1)) #TB
#    bPlot2DArrayOld = np.reshape(bPlotArray,(5,-1)) #r3
#    vStacked = np.dstack((1.-qPlotArray, wPlotArray))[0,:,:]
#    print vStacked

    pStartArray = np.linspace(0.1,1,num=1000)
    wStartArray = np.linspace(0,20,num=1000)

    grid_p, grid_w = np.meshgrid(pStartArray,wStartArray)

#    grid_q, grid_w = np.mgrid[0:1:1000j,0:20:1000j]
#    print grid_w
#    print grid_q

    bPlotArrayMod = bPlotArray.copy()
    bPlotArrayMod[bPlotArrayMod < -0.1] = 0.953
    print bPlotArrayMod

    bPlot2DArray = scipy.interpolate.griddata((1.-qPlotArray, wPlotArray), bPlotArrayMod, (grid_p, grid_w), method='cubic')
#    print bPlot2DArray
#yerr=errorBarArray
#    sc = plt.scatter(1.-qPlotArray, wPlotArray, c=bPlotArray, vmin=0, vmax=1, s=1296, cmap=cm, clip_on=False)
    im = plt.imshow(bPlot2DArray, vmin=0, vmax=1., interpolation='nearest',origin='lower',aspect='auto',extent=[0,plotWidth,0,plotHeight], cmap=cm)
    cbar = plt.colorbar(im)

    cbar.outline.set_linewidth(3.0)
    cbar.ax.set_ylabel(r'$b$',size=24,rotation=0,labelpad=19,verticalalignment="center")
    cbar.ax.tick_params(direction="inout", length=7, width=2, colors='k', top='off', right='on', labelsize=16)
#    cbar.ax.set_ylabel_size(20)





#    plt.plot(1.-qPlotArray, wPlotArray, '.', color='0.5', markersize=10, clip_on=False)
    plt.scatter(1.-qPlotArray, wPlotArray, facecolors='w', edgecolors='k', s=50, clip_on=False,zorder=100)

#    figNum += 1
#    fig = plt.figure(figNum,facecolor="white")
#    ax = plt.subplot()
#    im2 = plt.imshow(np.fliplr(bPlot2DArrayOld), vmin=bPlot2DArray.min(), vmax=bPlot2DArray.max(), interpolation='nearest',origin='lower',extent=[0,1,0,22.5],aspect='auto')
#    plt.colorbar(im)
#    plt.plot(1.-qPlotArray, wPlotArray, 'k.', markersize=50)

#vmin=0.0, vmax=1.0

#    line = ax.errorbar(qPlotArray, bPlotArray, yerr=errorBarArray, marker=markerList[typeIndex], markersize=markerSizeList[typeIndex], label="N={0}".format(NvalueCurrent), clip_on=False, linewidth=1.0, color=colourList[typeIndex])#linewidth=3.0, ls=''
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



#plt.ylim(0.,np.max(wPlotArray)*1.1)
#plt.xlim(0.,1.01)
plt.ylim(0.,plotHeight)
plt.xlim(0.,plotWidth)
#plt.xlabel(r'Probability of Occupation, $p$', fontsize=16)
#plt.ylabel(r'Disorder Strength, $w$', fontsize=16)
plt.xlabel(r'$p$', fontsize=24, labelpad=-3.)
plt.ylabel(r'$w$', fontsize=24)
#plt.title("Local Level Density", fontsize=18)
#plt.colorbar(sc)
#plt.colorbar(im)
[i.set_linewidth(3.0) for i in ax.spines.itervalues()]
for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(8.)
    tick.label1 = tick._get_text1()
for tick in ax.get_yaxis().get_major_ticks():
    tick.set_pad(10.)
    tick.label1 = tick._get_text1()
for label in ax.yaxis.get_ticklabels():
    label.set_verticalalignment('center')

ax.tick_params(direction="inout", length=10, width=2, colors='k', top='off', right='off', labelsize=20)

plt.tight_layout()

figureFilename = "/Users/jtcantin/Documents/phaseDiagram_N30_alphaINF.eps"
fig.savefig(figureFilename, format='eps', dpi=1200)


plt.show()




