#To run: python eigValCompTest.py numCompar LapackEvalFilenames... pythEigvalFilenames...

import sys
import numpy as np

numComparisons = int(sys.argv[1])

for i in range(2,numComparisons+2):
    LapackEvalFilename = sys.argv[i]
    pythEigvalFilename = sys.argv[i+numComparisons]


    pythEigvalFile = open(pythEigvalFilename, 'r')
    LapackEvalFile = open(LapackEvalFilename, 'r')

    pythEigval_Array = np.array(map(float, pythEigvalFile.readlines()))
    #print pythEigval_Array
    #print pythEigval_Array.shape

    LapackEval_Array = np.array(map(float, LapackEvalFile.readlines()))
    #print LapackEval_Array
    #print LapackEval_Array.shape

    absDiff = np.abs(pythEigval_Array - LapackEval_Array)
    relAbsDiff = absDiff / np.abs(pythEigval_Array)
    maxRelAbsDiff = np.max(relAbsDiff)

    LapackEvalFilename_short = LapackEvalFilename.split("/")[-1]
    pythEigvalFilename_short = pythEigvalFilename.split("/")[-1]
    
    print "------------------------"
    print "For files {0} and {1}".format(LapackEvalFilename_short, pythEigvalFilename_short)
    print "maxRelAbsDiff: ", maxRelAbsDiff

    pythEigvalFile.close()
    LapackEvalFile.close()


