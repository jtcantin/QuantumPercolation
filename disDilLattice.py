#
#  disDilLattice.py
#
#
# Created by Joshua Tyler Cantin on 2015-05-28.
#
#To run: python disDilLattice.py N occupationArrayFilename onsiteEnergyArrayFilename

import numpy as np
import sys

np.set_printoptions(threshold=10000,linewidth=2000,precision=6,suppress=False)

#Define a function to create conversion arrays from index to coordinates
def convArrays(dim,lengths):
    if (dim == 3):
        if len(lengths) != 3:
            sys.exit("ERROR: Need 3 lengths for 3D")

        Nx = lengths[0]
        Ny = lengths[1]
        Nz = lengths[2]

        indToCoorArray = []
        coorToInd3DArray = np.zeros((Nx,Ny,Nz))
        
        ind = 0
        
        for k in range(0,Nz):
            for j in range(0,Ny):
                for i in range(0,Nx):
                    indToCoorArray.append((i,j,k))
                    coorToInd3DArray[i,j,k] = ind

                    ind += 1

#        print coorToInd3DArray
#        print indToCoorArray
#        print Nx, Ny, Nz

        return coorToInd3DArray, indToCoorArray

    else:
        sys.exit("ERROR: Dim = %d not implemented" % dim)

#Parameters
t = 1.
alpha = 3
longRange = True

#Read in files
N = int(sys.argv[1])
occupationArrayFilename = sys.argv[2]
onsiteEnergyArrayFilename = sys.argv[3]

occupationArrayFile = open(occupationArrayFilename, 'r')

tester = 0
for line in occupationArrayFile:
    if (line[0] == "#") or (line.strip() == ""):
        continue
    
    if tester > 0:
        sys.exit("ERROR: More than one uncommented line encountered.")

    tmpArray = map(int, line.split())
    nSitesOc, occupationArray = tmpArray[0], tmpArray[1:]

#    print nSitesOc
#    print occupationArray

    tester += 1

occupationArrayFile.close()

onsiteEnergyArrayFile = open(onsiteEnergyArrayFilename, 'r')

print " "
print " "

tester = 0
for line in onsiteEnergyArrayFile:
    if (line[0] == "#") or (line.strip() == ""):
        continue
    
    if tester > 0:
        sys.exit("ERROR: More than one uncommented line encountered.")

    nSitesEng = int(line.split()[0])
    onsiteEnergyArray = map(float, line.split()[1:])


#    print nSitesEng
#    print onsiteEnergyArray

    tester += 1

onsiteEnergyArrayFile.close()

if nSitesOc != nSitesEng:
    sys.exit("ERROR: Number of sites in the two files do not match.")

#Build Conversion arrays
print "Side Length: %d" % N
print "Lattice Size: %d" % N**3
print "Number of Occupied Sites: %d" % nSitesOc

coorToInd3DArray, indToCoorArray = convArrays(3,(N,N,N))

print "Building Hamiltonian."

#Build Hamiltonian
hami2DArray = np.zeros((nSitesOc,nSitesOc))

for i in range(0, nSitesOc):
    for j in range(0, nSitesOc):
        if i == j:
            hami2DArray[i,j] = onsiteEnergyArray[i]
        
        else:
            #Get site coordinates
            iVec = np.array(indToCoorArray[occupationArray[i]])
            jVec = np.array(indToCoorArray[occupationArray[j]])
            
            #Get distance between sites
            diffVec = iVec - jVec
            dist = np.linalg.norm(diffVec)
#            if i==13:
#                print i
#                print j
#                print iVec
#                print jVec
#                print dist

            if longRange:
                #Assign element
                hami2DArray[i,j] = t/(dist**alpha)
            
            else: #Then Tight Binding model
                if (dist < 1.001): #True if nearest neighbour
                    hami2DArray[i,j] = t
#                    print "triggered"

print "Hamiltonian built."
print "Diagonalizing..."
#np.savetxt("test.csv",hami2DArray,delimiter=",")
#print hami2DArray
#print np.max(np.abs(hami2DArray.T - hami2DArray))
#print indToCoorArray

#Diagonalize the Hamiltonian
#eval, evec = np.linalg.eig(hami2DArray)
eval = np.linalg.eigvalsh(hami2DArray)

eval_srtd = np.sort(eval)

#sort the eigenvectors, so that they are in order of increasing energy
#evec_srtd = np.copy(evec[:,eval.argsort()])

print "Eigenvalues calculated and sorted."

#print " "
#print "Eigenvalues:"
#for eval_item in eval_srtd:
#    print eval_item

eigvalFilename = occupationArrayFilename[:-4]+"_eigvals.txt"
np.savetxt(eigvalFilename,eval_srtd,delimiter=',')
print "Eigenvalues saved to :", eigvalFilename










print " "





























