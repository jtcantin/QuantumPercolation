#This is to test whether numpy.random.permuation rearranges the
# coefficients as needed.
import sys
import numpy as np
import matplotlib.pyplot as plt

N = int(sys.argv[1])
M = int(sys.argv[2])

initArray = np.arange(N)+1
countArray = []
sumArray = np.zeros(N)
permedArrayOld = np.zeros(N)
corrArray = np.zeros(N)

for i in range(0,M):
    permedArray = np.random.permutation(initArray)
    corrArray = corrArray + permedArrayOld * permedArray
#    print permedArrayOld
#    print permedArray
#    print corrArray
    countArray.append(permedArray[0])
    sumArray = sumArray + permedArray
    permedArrayOld = permedArray.copy()

sumArray = sumArray / M
corrArray = (corrArray/M) - ((sumArray)**2)
#print corrArray

#print countArray
print "Expected Avg: ", (N+1)/2.
print "Max Avg: ", np.max(sumArray)
print "Mean Avg: ", np.mean(sumArray)
print "Min Avg: ", np.min(sumArray)

#Plot
figNum = 0

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
datas = ax.hist(countArray, np.arange(N+1)+1./2,normed=True)
datas2, = ax.plot(np.arange(N),(1./N*np.ones(N)),color='red')

plt.ylim(0,2./N)
#plt.xlim(0,N)
plt.xlabel("Number", fontsize=16)
plt.ylabel("Probability", fontsize=16)
plt.title("Number Histogram", fontsize=18)

figNum += 1
fig = plt.figure(figNum)
ax = plt.subplot()
datas2, = ax.plot(initArray, corrArray)

#plt.ylim(0.9,1.1)
#plt.xlim(0,N)
plt.xlabel("Number", fontsize=16)
plt.ylabel("Relative Correlation", fontsize=16)
plt.title("Two Step Correlation", fontsize=18)

plt.show()