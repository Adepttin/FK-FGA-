import subprocess
import numpy as np
from shutil import copyfile
from routines import *

'''
Find chemical potential mu to fix wanted occupation of mobile electrons
'''

#inverste temperature
beta = 10.
#interaction strength
U = 0.5
#occupation of localized f electrons
p1 = 0.5
#DF lattice size
nk = 6
#DMFT lattice size
nkDMFT = 128
#number of DF Matsubara frequencies
nv = 20
#number of DMFT Matsubara frequencies - take care to use enough!
nvDMFT = 100
#number of iterations in DMFT
DMFTiter = 40

#wanted occupation of mobile electrons
n = 0.15

#search interval for mu
mu_a = 0.
mu_b = -1.

#call findMu
mu, sigma, occupation = findMu(mu_a,mu_b,n,1e-5,beta,U,p1,nkDMFT,nk,nvDMFT,nv,DMFTiter)

#print found mu and corresponding sigma (deviation from half filling) and occupation
print("---------------------")
print("Results:")
print(mu, sigma, occupation)
print("---------------------")

	
	
	
	
	
	
	

