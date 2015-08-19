#Written by Tianrui Xu
import numpy as np
import os

def rand_vac(size, ptg, num_dis, out_dir):
    out_file='rd3d'+str(int(size))+str(int(ptg))
    pct=ptg*0.01
    ttl=size*size*size
    N=int(pct*ttl)
    print N
    rng=np.arange(0,ttl)
    x=rng
    os.system('mkdir -p '+out_dir)
    for t in range(0,num_dis):
        a=np.random.permutation(rng)[:N]
        a=np.sort(a)
        f = open(out_dir+'/'+out_file+str(int(t))+'.txt','w')
        f.write(str(N))
        f.write('\t')
        for i in range(0,len(a)):
            f.write(str(a[i]))
            f.write('\t')
        f.close()

def rand_occ(size, ptg, num_dis, out_dir):
    out_file='rd3do'+str(int(size))+str(int(ptg))
    #pct=ptg*0.01       #<-- leave this for now
    #ttl=size*size*size #<-- leave this for now
    #N=int((1-pct)*ttl) #<-- leave this for now
    #N=ttl-int(pct*ttl) #<-- leave this for now
    ptg=100-ptg
    pct=ptg*0.01
    ttl=size*size*size
    N=int(pct*ttl)#<--- well this line would force the function returns an integer number of numbers, would not be too much of a trouble for most of the time
    print "number of total lattice sites = "+str(ttl)
    print "number of sites occupied = "+str(N)
    rng=np.arange(0,ttl)
    x=rng
    os.system('mkdir -p '+out_dir)
    for t in range(0,num_dis):
        a=np.random.permutation(rng)[:N]
        a=np.sort(a)
        f = open(out_dir+'/'+out_file+str(int(t))+'.txt','w')
        f.write(str(N))
        f.write('\t')
        for i in range(0,len(a)):
            f.write(str(a[i]))
            f.write('\t')
        f.close()

def rand_onsite(size, ptg, minw, maxw, num_dis, out_dir):
    out_file='rd3dw'+str(int(size))+str(int(ptg))
    #pct=ptg*0.01       #<-- leave this for now
    #ttl=size*size*size #<-- leave this for now
    #N=int((1-pct)*ttl) #<-- leave this for now
    #N=ttl-int(pct*ttl) #<-- leave this for now
    ptg=100-ptg
    pct=ptg*0.01
    ttl=size*size*size
    N=int(pct*ttl)#<--- well this line would force the function returns an integer number of numbers, would not be too much of a trouble for most of the time
    print "number of total lattice sites = "+str(ttl)
    print "number of sites occupied = "+str(N)
    os.system('mkdir -p '+out_dir)
    for t in range(0,num_dis):
        w = np.random.uniform(minw,maxw,N)
        print 'infimum = '+str(w.min())+' supremum = '+str(w.max())
        f = open(out_dir+'/'+out_file+str(int(t))+'.txt','w')
        f.write(str(N))
        f.write('\t')
        for i in range(0,len(w)):
            f.write(str(w[i]))
            f.write('\t')
        f.close()
