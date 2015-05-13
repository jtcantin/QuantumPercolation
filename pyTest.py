import numpy as np
import sys

np.set_printoptions(threshold=10000,linewidth=2000,precision=4,suppress=False)

N = int(sys.argv[1])
E = float(sys.argv[2])
t = float(sys.argv[3])

Hami = np.zeros((N,N))

for i in range(0,N):
    for j in range(0,N):
        if (i==j):
            Hami[i][j] = E

        elif (np.abs(i-j) == 1):
            Hami[i][j] = t

eval_Array, evec_2DArray = np.linalg.eig(Hami)

eval_Array_srtd = np.sort(eval_Array)
    
#sort the eigenvectors, so that they are in order of increasing energy
evec_2DArray_srtd = np.copy(evec_2DArray[:,eval_Array.argsort()])

print "Eigenvalues: "
#print eval_Array_srtd
for i in range(0,N):
    print "%.30f" % eval_Array_srtd[i]