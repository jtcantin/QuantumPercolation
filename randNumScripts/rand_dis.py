#Written by Tianrui Xu
import rand_dis_fcn as rdf

#All arguments except 'dir' are numbers

dir='/Users/tianruixu/Documents/code/disorders/disorders' #<-- directory to store random files, needs to be a string

size=51 # <-- size of the lattice per dimension, integer
ptg_vac=90 #<-- percentage of vacancy, prefer integer, otherwise might need to modify 'out_file' in functions
num_dis=1 #<-- number of disorders generated, integer
wmin=-1 #<-- infimum of on-site energy
wmax=1 #<-- supremum of on-site energy


#Although the functions take vacancy percentage, they will return occupation stuff. The reason for this is that I used random vacancy more often when I was doing my previous project and I just need to add one more line to my previous script to return random occupied sites... so the functions still take vacancy percentage as its argument... :P

rdf.rand_occ(size,ptg_vac,num_dis,dir) #<-- random occupied site number
rdf.rand_onsite(size, ptg_vac, wmin, wmax, num_dis, dir) #<-- random on-site energy
