#test stuff in t_in.txt
import os
import sys

print 'take arguments: '
print sys.argv

c_num=sys.argv[1]
size=sys.argv[2]
pvac=sys.argv[3]

os.system('cp /home/tianruix/QuantumPercolation/TestCode/t_in.txt ./dis_in.txt')
os.system('cp /home/tianruix/QuantumPercolation/TestCode/rd3d?'+str(size)+str(pvac)+'0.txt .')
chSize='sed -i "s/TrSize/'+str(size)+'/g" dis_in.txt'
os.system(chSize)
chPvac='sed -i "s/TrPvac/'+str(pvac)+'/g" dis_in.txt'
os.system(chPvac)

os.system('date')
os.system('./QmPerc_LINUX 0')
os.system('date')

