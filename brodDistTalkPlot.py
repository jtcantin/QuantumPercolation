import numpy as np
import sys
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import distFit

SvalsArray = np.linspace(0,7,1000)
bDistWignerArray = distFit.brodyDist(SvalsArray, 1.0)
bDistPoissonArray = distFit.brodyDist(SvalsArray, 0.0)
bDistMiddleArray = distFit.brodyDist(SvalsArray, 0.5)

figNum = 0

figNum += 1
fig = plt.figure(figNum,facecolor="white")
ax = plt.subplot()
#    binEdgesArrayArray_qN_srtd[distNum][:-1]
line, = ax.plot(SvalsArray, bDistPoissonArray,'b-',linewidth=3.0, clip_on=False, label="b = 0.0, Poissonian")
line2, = ax.plot(SvalsArray, bDistWignerArray,'r-',linewidth=3.0, clip_on=False, label="b = 1.0, Wigner")
line3, = ax.plot(SvalsArray, bDistMiddleArray,'k-',linewidth=3.0, clip_on=False, label="b = 0.5")
#line, = ax.plot(NArray_srtd, bValArray_srtd,'k-',marker=".",markersize=20,linewidth=2.0)#, label="Local Level Density")
#line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
ax.legend()

# Remove plot frame
#ax.set_frame_on(False)



plt.ylim(0,1.0)
#plt.xlim(0,N)
plt.xlabel("Energy Level Spacing, S", fontsize=16)
plt.ylabel("P(S)", fontsize=16)
#    plt.title("DistNum: {0}, bVal: {1:4f}".format(distNum,data2DArray_qN_srtd[3][distNum]), fontsize=18)

plt.show()



