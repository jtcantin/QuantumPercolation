import sys
import numpy as np

nOccSites = int(sys.argv[1])
hamiArrayFilename = sys.argv[2]

hamiArrayFile = open(hamiArrayFilename, 'r')

hamiArray = map(float, hamiArrayFile.readline().split())

hami2DArray = np.resize(hamiArray, (nOccSites,nOccSites))

#print hami2DArray

np.savetxt("trxTest.csv",hami2DArray,delimiter=",")