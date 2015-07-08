import numpy as np
import scipy.special as spc
import sys
import matplotlib.pyplot as plt


#########################
#########################
#Functions
#########################
#########################
def brodyDist(S,q):
    q1 = 1.+q
    
    beta = np.power(spc.gamma((2.+q)/q1), q1)
    
    alpha = q1*beta

    factor1 = alpha*np.power(S, q)

    exponent = -1.*beta*np.power(S,q1)

    factor2 = np.exp(exponent)

    return factor1*factor2

#########################
#########################
#Program
#########################
#########################

numPoints = int(sys.argv[1])
b = float(sys.argv[2])
bNum = float(sys.argv[3])

S = np.linspace(0, 10, numPoints)

bArray = np.linspace(0,1,bNum)

brodyArray = brodyDist(S,b)

brodyBList = []

for i in range(0, bArray.size):
    brodyBList.append(brodyDist(S,bArray[i]))

#########################
#########################
#Plot
#########################
#########################
figNum = 0

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
for i in range(0, bArray.size):
    labelString = "b = %f" % bArray[i]
    line, = ax.plot(S, brodyBList[i], label=labelString)


#line, = ax.plot(S, brodyArray)#, label="Local Level Density")
#line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
ax.legend()

plt.ylim(0,1.0)
#plt.xlim(0,N)
plt.xlabel(r'$S$', fontsize=16)
plt.ylabel(r'$P(S)$', fontsize=16)
plt.title("Brody Distribution", fontsize=18)

#########################################
figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
#for i in range(0, bArray.size):
#    labelString = "b = %f" % bArray[i]
#    line, = ax.plot(S, brodyBList[i], label=labelString)

labelString = "b = %f" % b
line, = ax.plot(S, brodyArray, label=labelString)
#line2, = ax.plot(eigvals, rhoGauss, color='red', marker="o", label="Gaussian Broadening Method")
ax.legend()

plt.ylim(0,1.0)
#plt.xlim(0,N)
plt.xlabel(r'$S$', fontsize=16)
plt.ylabel(r'$P(S)$', fontsize=16)
plt.title("Brody Distribution", fontsize=18)

plt.show()
