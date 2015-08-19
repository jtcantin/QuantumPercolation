import rand_dis_fcn as rdf
import numpy as np
import time

N = 100
q100 = 0
num_dis = 5
w = 72.3546
dir = "randomTestFiles"
#print time.time()
seedVal = int(100*time.time()) % 4294967295


wmin = -w/2.
wmax = w/2.

#Set seed
np.random.seed(seedVal)

#Use Tianrui's code to generate disorder
rdf.rand_occ(N,q100,num_dis,dir) #<-- random occupied site number
rdf.rand_onsite(N, q100, wmin, wmax, num_dis, dir) #<-- random on-site energy

#Read in Tianrui's Disorders
occTRXlist = []
onSiteTRXlist = []

for dis in range(0,num_dis):
    occRandFile = open(dir + "/rd3do{0}{1}{2}.txt".format(N,q100,dis),'r')
    onSiteRandFile = open(dir + "/rd3dw{0}{1}{2}.txt".format(N,q100,dis),'r')
    
    occRandLine = occRandFile.readline()
    occRandLineSplit = occRandLine.split()
#    numOccSites = int(occRandLineSplit[0])
    occRandLineSplitStripped = map(lambda s: s.strip(),occRandLineSplit[1:])
    occRandArray = np.array(map(int,occRandLineSplitStripped))
    occTRXlist.append(occRandArray)

    onSiteRandLine = onSiteRandFile.readline()
    onSiteRandLineSplit = onSiteRandLine.split()
    #    numOnSiteSites = int(onSiteRandLineSplit[0])
    onSiteRandLineSplitStripped = map(lambda s: s.strip(),onSiteRandLineSplit[1:])
    onSiteRandArray = np.array(map(float,onSiteRandLineSplitStripped))
    onSiteTRXlist.append(onSiteRandArray)

    occRandFile.close()
    onSiteRandFile.close()

#Generate my own disorders, using the same numpy routines as Tianrui:

#Reset Seed
np.random.seed(seedVal)

#Prelim work
numLatticeSites = N*N*N
numSites = int(numLatticeSites*(1.-(q100/100.)))

#Generate occupation disorder first
occJTClist = []

for dis in range(0,num_dis):
    randArray = np.random.permutation(numLatticeSites)[:numSites]
    occJTClist.append(np.sort(randArray))

#Generate on site disorder next
onSiteJTClist = []

for dis in range(0,num_dis):
    randArray = np.random.uniform(wmin,wmax,numSites)
    onSiteJTClist.append(randArray)

#Compare results
occMaxAbsDiffList = []
onSiteMaxAbsDiffList = []
onSiteMaxRelDiffList = []

for dis in range(0,num_dis):
    #Occupation disorder
    absDiff = np.abs(occJTClist[dis] - occTRXlist[dis])
    maxAbsDiff = np.max(absDiff)
    print "Max abs diff for occupied sites, dis{0}: {1}".format(dis, maxAbsDiff)
    occMaxAbsDiffList.append(maxAbsDiff)

    #On site disorder
    absDiff = np.abs(onSiteJTClist[dis] - onSiteTRXlist[dis])
    maxAbsDiff = np.max(absDiff)
    
    relDiff = np.abs(absDiff/onSiteJTClist[dis])
    maxRelDiff = np.max(relDiff)
    
    print "Max abs diff for on site disorder, dis{0}: {1}".format(dis, maxAbsDiff)
    onSiteMaxAbsDiffList.append(maxAbsDiff)
    
    print "Max rel diff for on site disorder, dis{0}: {1}".format(dis, maxRelDiff)
    onSiteMaxRelDiffList.append(maxRelDiff)

print "----"
print "Max abs diff for occupied sites: ", np.max(np.array(occMaxAbsDiffList))
print "Max abs diff for on site disorder: ", np.max(np.array(onSiteMaxAbsDiffList))
print "Max rel diff for on site disorder: ", np.max(np.array(onSiteMaxRelDiffList))















