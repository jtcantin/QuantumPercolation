#Written by Tianrui Xu
#Modified by Joshua Cantin on 2015/06/08
import rand_dis_fcn as rdf
import datetime

#All arguments except 'dir' are numbers

dir='rndDatafiles/' #<-- directory to store random files, needs to be a string


size=10 # <-- size of the lattice per dimension, integer
ptg_vac=01 #<-- percentage of vacancy, prefer integer, otherwise might need to modify 'out_file' in functions
num_dis=1 #<-- number of disorders generated, integer
wmin=-0.1 #<-- infimum of on-site energy
wmax=0.1 #<-- supremum of on-site energy

dirAdd = ("n%03dq%02dw%03d" % (size,ptg_vac,wmax*2))
#print dirAdd
dir = dir + "files_" + dirAdd + "_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") #Add current date and time to differentiate folders

#print dir

#Although the functions take vacancy percentage, they will return occupation stuff. The reason for this is that I used random vacancy more often when I was doing my previous project and I just need to add one more line to my previous script to return random occupied sites... so the functions still take vacancy percentage as its argument... :P

rdf.rand_occ(size,ptg_vac,num_dis,dir) #<-- random occupied site number
rdf.rand_onsite(size, ptg_vac, wmin, wmax, num_dis, dir) #<-- random on-site energy
