#To run: python eigValCompTest.py numCompar LapackEvecFilenames... pythEigvecFilenames...

import sys
import numpy as np

np.set_printoptions(threshold=10000,linewidth=2000,precision=6,suppress=False)

numComparisons = int(sys.argv[1])

for i in range(2,numComparisons+2):
    LapackEvecFilename = sys.argv[i]
    pythEigvecFilename = sys.argv[i+numComparisons]


    pythEigvecFile = open(pythEigvecFilename, 'r')
    LapackEvecFile = open(LapackEvecFilename, 'r')

    pythEigvec_List = []
    
    #Skip first three comment lines
    pythEigvecFile.readline()
    pythEigvecFile.readline()
    pythEigvecFile.readline()
    
    for line in pythEigvecFile:
        splitLine = line.split()
        strip_list = map(lambda it: it.strip(), splitLine)
        pythEigvec_List.append(np.array(map(float, strip_list[1:])))
#        print pythEigvec_List
#        raw_input()

    pythEigvec_Array = np.vstack(pythEigvec_List)
    #print pythEigvec_Array
    print pythEigvec_Array.shape

    LapackEvec_List = []

    #Skip first two comment lines
    LapackEvecFile.readline()
    LapackEvecFile.readline()
    
    for line in LapackEvecFile:
        splitLine = line.split()
        strip_list = map(lambda it: it.strip(), splitLine)
        LapackEvec_List.append(np.array(map(float, strip_list[1:])))

    LapackEvec_Array = np.vstack(LapackEvec_List)
    #print LapackEvec_Array
    print LapackEvec_Array.shape

    #Make sure all vectors have the same global phase
    pythPhase = pythEigvec_Array[0,:] > 0
    lapackPhase = LapackEvec_Array[0,:] > 0

    phasesToChange = np.logical_xor(pythPhase,lapackPhase)
    multArray = np.ones(pythEigvec_Array[0].size) - 2*phasesToChange

#    print pythPhase[0:10]
#    print lapackPhase[0:10]
#    print phasesToChange[0:10]
#    print multArray[0:10]

#    print pythEigvec_Array[0:10,0:10]

#    old = np.copy(pythEigvec_Array[0:10,0:10])

    pythEigvec_Array = pythEigvec_Array[:] * multArray

#    print pythEigvec_Array[0:10,0:10]

#    print pythEigvec_Array[0:10,0:10] + old

    absDiff = np.abs(pythEigvec_Array - LapackEvec_Array)
    relAbsDiff = absDiff / np.abs(pythEigvec_Array)
    maxRelAbsDiff = np.max(relAbsDiff)

    LapackEvecFilename_short = LapackEvecFilename.split("/")[-1]
    pythEigvecFilename_short = pythEigvecFilename.split("/")[-1]
    
    print "------------------------"
    print "For files {0} and {1}".format(LapackEvecFilename_short, pythEigvecFilename_short)
    print "maxRelAbsDiff: ", maxRelAbsDiff
#    print np.min(relAbsDiff)

    pythEigvecFile.close()
    LapackEvecFile.close()


