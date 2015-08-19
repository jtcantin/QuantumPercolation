import sys
import os
import datetime
import numpy as np
import rand_dis_fcn as rdf


runNum = 5
numThreads = 8

alphaList = [3]
qList = [90]
wList = [20]
Nlist = [103]
numDisorders = 1

#q=25,40 transtion line: wList = [0,2,4,5,6,7,8,8.5,9,9.5,10,10.5,11,12,13,15,20]
#q=0 transition line (tentative): [1,5,8,10,12,13,14,15,15.5,16,16.5,17,18,19,20,22,24,28]
#w=4 transition line: [90, 80, 70, 60, 55, 50, 45, 40, 35, 30, 20, 10, 0]
#Initial phase diagram spread: N=30, qList = [90, 75, 50, 25,0], wList = [0, 5, 10, 15, 20], numDisorders = 3

def dis_in_make(N,rng,q,filename):
    
    if np.isinf(rng):
        numNN = 1
        gamma = 3
    else:
        numNN = "Long"
        gamma = rng
    
    dis_inFile = open(filename,'w')

    dis_inFile.write("!Size -- has to be an integer, the actual size of the system is size^3\n")
    dis_inFile.write("{0}\n".format(N))
    
    dis_inFile.write("!Range -- xNN/ 'Long'\n")
    dis_inFile.write("{0}\n".format(numNN))

    dis_inFile.write("!Gamma -- tunnelling amplitude ~d^(-Gamma), double\n")
    dis_inFile.write("{0}\n".format(gamma))

    dis_inFile.write("!Z-hop -- 0: t_Z=t; 1: t_Z=2t\n")
    dis_inFile.write("0\n")

    dis_inFile.write("!Percent of Vacancy -- percentage, i.e. x%\n")
    dis_inFile.write("{0}\n".format(q))

    dis_inFile.write("!end of the input file\n")
    dis_inFile.write("\n")

    dis_inFile.close()

bigRunDir = "run{0}".format(runNum)

bashFilename = "%s/runCode%d.sh" % (bigRunDir, runNum)

bashFile = open(bashFilename, 'w')

tmpDir = "/tmp/jtc/QPerc/run%d/" % runNum

#BASH: Move relevant files to tmp directory
bashFile.write("cp -r * %s \n" % tmpDir)
bashFile.write("cd %s \n" % tmpDir)
bashFile.write("\n")

#BASH: Verify current directory
bashFile.write("echo $PWD \n")
bashFile.write("\n")

#BASH: Set and display number of threads
bashFile.write("export OMP_NUM_THREADS=%d \n" % numThreads)
bashFile.write("export MKL_NUM_THREADS=%d \n" % numThreads)
bashFile.write("\n")

bashFile.write("echo OMP_NUM_THREADS \n")
bashFile.write("echo $OMP_NUM_THREADS \n")
bashFile.write("\n")
bashFile.write("echo MKL_NUM_THREADS \n")
bashFile.write("echo $MKL_NUM_THREADS \n")
bashFile.write("\n")

#Execution Commands
#N = 30
#q = 60
#w = 13
#alpha = 3.2

for alpha in alphaList:
    for q in qList:
        for w in wList:
            for N in Nlist:

                if np.isinf(alpha):
                    rng = "TB"
                else:
                    integerPart = int(alpha)
                    decimalPart = alpha - integerPart
                    rng = "{0}-{1}".format(integerPart, str(decimalPart)[2:])

                #Make disorder files and dataDir
                currentTimeString = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                if np.abs((int(w) - w)) < 1E-16:
                    dataDirStem = "data_N{0:03d}_q{1:02d}_w{2:03d}_r{3}_{4}".format(N,q,w,rng,currentTimeString)
                else:
                    dataDirStem = "data_N{0:03d}_q{1:02d}_w{2:03.2f}_r{3}_{4}".format(N,q,w,rng,currentTimeString)
                        
                        
                dataDir = "{0}/{1}".format(bigRunDir,dataDirStem)

                wmin = -w/2.
                wmax = w/2.

                rdf.rand_occ(N,q,numDisorders,dataDir) #<-- random occupied site number
                rdf.rand_onsite(N, q, wmin, wmax, numDisorders, dataDir) #<-- random on-site energy

                #Make input file
                inputFile = "dis_in_N{0}_q{1}_w{2}_r{3}.txt".format(N, q, w, rng)
                inputFileLoc = "{4}/dis_in_N{0}_q{1}_w{2}_r{3}.txt".format(N, q, w, rng, dataDir)
                #inputFile = "dis_in_N{0}_q{1}.txt".format(N, q)

                dis_in_make(N,alpha,q,inputFileLoc)

                #Copy executable files to dataDir
                os.system("cp QmPerc_LINUX_TRX {0}/".format(dataDir))

                #BASH: cd into dataDirStem
                bashFile.write("cd %s \n" % dataDirStem)
                
                bashFile.write("echo \"-----------------Starting calculations for {0}---------------------\"\n".format(dataDirStem))
                bashFile.write("\n")

                #BASH: Run calculations
                for disNum in range(0,numDisorders):
                    bashFile.write("echo \"----------------\" \n")
                    bashFile.write("echo \"N%d q%d w%d r%s Dis%d\" \n" % (N, q, w, rng, disNum) )
                    bashFile.write("/usr/bin/time ./QmPerc_LINUX_TRX {0} {1} \n".format(inputFile, disNum))
                    bashFile.write("\n")

                bashFile.write("echo \"-----------------Calculations for {0} completed---------------------\"\n".format(dataDirStem))
                bashFile.write("\n")

                #BASH: cd back out of dataDir
                bashFile.write("cd .. \n")
                bashFile.write("\n")

##################################

bashFile.write("echo \"-----------------All calculations completed---------------------\"\n")
bashFile.write("\n")

#BASH: Copy everything to data root directory
dataRootDir = "/home/jtcantin/code/QuantumPercolation/srcTRX/data/"
#dataStorageDir = dataRootDir + dataDir


bashFile.write("cp -r /tmp/jtc/QPerc/run{0}/* {1} \n".format(runNum, dataRootDir))
bashFile.write("\n")

bashFile.write("mv {0}runCode{1}.sh {0}runCode{1}_{2}.sh \n".format(dataRootDir, runNum, currentTimeString))
bashFile.write("\n")

bashFile.write("echo \"-------------------Copying completed-----------------------\"\n")
bashFile.write("\n")
         
#BASH: Email saying the job is done
bashFile.write("echo \"Run{0} has finished.\" | mail -s \"Calculation Finished\" jtcantin@chem.ubc.ca \n".format(runNum))

bashFile.close()

#Make the bash script executable
os.system("chmod 744 {0}".format(bashFilename))

print "Files created."
print "Execution command:"
print "./runCode{0}.sh &> {1}log{2}.txt".format(runNum,tmpDir,currentTimeString)


