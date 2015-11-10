# To run: name.py [inputDir]


import numpy as np
import sys
from os import listdir
from os.path import isfile, join
#import matplotlib
#matplotlib.use('Qt4Agg')
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

data2DArray_q_srtd_indices = np.argsort(data2DArray[1])

data2DArray_q_srtd = np.copy(data2DArray[:,data2DArray_q_srtd_indices])

#print data2DArray_q_srtd

first_qIndices = np.argsort(data2DArray_q_srtd[0,(np.abs(data2DArray_q_srtd[1]-data2DArray_q_srtd[1][0]) < 1E-16)])
second_qIndices = np.argsort(data2DArray_q_srtd[0,(np.abs(data2DArray_q_srtd[1]-data2DArray_q_srtd[1][-1]) < 1E-16)])+np.max(first_qIndices)+1

#print first_qIndices
#print second_qIndices

data2DArray_qN_srtd_indices= np.concatenate((first_qIndices,second_qIndices))

data2DArray_qN_srtd = np.copy(data2DArray_q_srtd[:,data2DArray_qN_srtd_indices])
print data2DArray_qN_srtd[:,0]
print data2DArray_qN_srtd.shape

#Sort distribution Arrays
distArrayArray_q_srtd = np.copy(distArrayArray[data2DArray_q_srtd_indices])
distArrayArray_qN_srtd = np.copy(distArrayArray_q_srtd[data2DArray_qN_srtd_indices])
binEdgesArrayArray_q_srtd = np.copy(binEdgesArrayArray[data2DArray_q_srtd_indices])
binEdgesArrayArray_qN_srtd = np.copy(binEdgesArrayArray_q_srtd[data2DArray_qN_srtd_indices])

#Generate Brody Distributions
bDistArrayList = []
SValsArrayList = []
for distNum in range(0, distArrayArray_qN_srtd.size):
    SValsArrayList.append(binEdgesArrayArray_qN_srtd[distNum][:-1] + (binEdgesArrayArray_qN_srtd[distNum][1]-binEdgesArrayArray_qN_srtd[distNum][0])/2.)

    bDistArrayList.append(distFit.brodyDist(SValsArrayList[distNum], data2DArray_qN_srtd[3][distNum]))


#data2DArray_qN_srtd = np.sort(data2DArray,axis=1,order=[1,0]


#Sort one array given the sorting of another array
#NArray_srtd = np.sort(NArray)
#bValArray_srtd = np.copy(bValArray[NArray.argsort()])
#
#print NArray_srtd
#print bValArray_srtd

#########################
#########################
#Plot
#########################
#########################
figNum = 0

#######
#b vs N
#######

first_q_bValArray = data2DArray_qN_srtd[3][:np.max(first_qIndices)+1]
second_q_bValArray = data2DArray_qN_srtd[3][np.max(first_qIndices)+1:]

first_q_NArray = data2DArray_qN_srtd[0][:np.max(first_qIndices)+1]
second_q_NArray = data2DArray_qN_srtd[0][np.max(first_qIndices)+1:]


figNum += 1
fig = plt.figure(figNum,facecolor="white")
ax = plt.subplot()
line, = ax.plot(first_q_NArray, first_q_bValArray,'k-',marker="D", markersize=10,linewidth=3.0, label="Diffusive Regime", clip_on=False)
line2, = ax.plot(second_q_NArray, second_q_bValArray,'k-',marker=".",markersize=20,linewidth=3.0, label="Localized Regime", clip_on=False)
#line, = ax.plot(NArray_srtd, bValArray_srtd,'k-',marker=".",markersize=20,linewidth=2.0)#, label="Local Level Density")
#line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
#ax.legend()

# Remove plot frame
#ax.set_frame_on(False)
#figTop = plt.get_current_fig_manager()

[i.set_linewidth(3.0) for i in ax.spines.itervalues()]
plt.ylim(0,1.)
#plt.xlim(0,N)
plt.xlabel(r'$N$', fontsize=24)
plt.ylabel(r'$b$', fontsize=24)
#plt.title("Local Level Density", fontsize=18)
fig.canvas.set_window_title('Size Convergence')
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
legend = ax.legend()
legend.get_frame().set_linewidth(2.0)
#fig.savefig('SizeConvergence.svg', format='svg', dpi=1200)
fig.savefig('SizeConvergence.eps', format='eps', dpi=1200)
#fig.savefig('SizeConvergence.svg', format='svg', dpi=1200)

#######
#Distributions
#######

first_q_bValArray = data2DArray_qN_srtd[3][:np.max(first_qIndices)+1]
second_q_bValArray = data2DArray_qN_srtd[3][np.max(first_qIndices)+1:]

first_q_NArray = data2DArray_qN_srtd[0][:np.max(first_qIndices)+1]
second_q_NArray = data2DArray_qN_srtd[0][np.max(first_qIndices)+1:]

for distNum in range(0, distArrayArray_qN_srtd.size):

    figNum += 1
    fig = plt.figure(figNum,facecolor="white")
    ax = plt.subplot()
#    binEdgesArrayArray_qN_srtd[distNum][:-1]
    line, = ax.plot(SValsArrayList[distNum], distArrayArray_qN_srtd[distNum],'ko', markersize=6, clip_on=False,label="Numerical Results")
    line2, = ax.plot(SValsArrayList[distNum], bDistArrayList[distNum],'r-',linewidth=2.0, clip_on=False, label="Brody Distribution Fit",zorder=100)
    #line, = ax.plot(NArray_srtd, bValArray_srtd,'k-',marker=".",markersize=20,linewidth=2.0)#, label="Local Level Density")
    #line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
    ax.legend()

    # Remove plot frame
    #ax.set_frame_on(False)


    [i.set_linewidth(3.0) for i in ax.spines.itervalues()]
    plt.ylim(0,1.0)
    #plt.xlim(0,N)
    plt.xlabel(r'$S$', fontsize=24)
    plt.ylabel(r'$P(S)$', fontsize=24)
#    plt.title("DistNum: {0}, bVal: {1:4f}".format(distNum,data2DArray_qN_srtd[3][distNum]), fontsize=18)
    fig.canvas.set_window_title("DistNum: {0}, bVal: {1:2f}, N:{2:d},q:{3},w:{4}".format(distNum,data2DArray_qN_srtd[3][distNum],int(data2DArray_qN_srtd[0][distNum]),data2DArray_qN_srtd[1][distNum],data2DArray_qN_srtd[2][distNum]))
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
    legend = ax.legend()
    legend.get_frame().set_linewidth(2.0)

    
    if int(data2DArray_qN_srtd[0][distNum]) == 35:
        figureFilename = "N{0:d}_q{1}_w{2}.eps".format(int(data2DArray_qN_srtd[0][distNum]),data2DArray_qN_srtd[1][distNum],data2DArray_qN_srtd[2][distNum])
        fig.savefig(figureFilename, format='eps', dpi=1200)

#figTop.canvas.manager.window.raise_()
plt.show()








