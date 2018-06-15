import subprocess
import numpy as np
from shutil import copyfile
from routines import *

'''
Find mu so that occupation of mobile electrons is n=0.15
'''

beta = 10.
#mu = -0.4
U = 0.5
Ef = 1.
p1 = 0.5
nk = 6
nkDMFT = 128
nv = 20
nvDMFT = 100 #200
DMFTiter = 40

dT = 0.005
Tstart = 0.1
NT = 21

#Temperatures from 0.1 to 0.02
Temperatures = [Tstart - dT*i for i in range(0,NT-1)]
#Temperatures.append(0.0001)

#print(Temperatures)

#calculating corresponding betas
Betas = [1./T for T in Temperatures]

#print(Betas)

muList = np.zeros((len(Betas)),dtype=float)

mu_a = 1.
mu_b = -2.

n=0.15


for i in range(0,len(Betas)):

	#for mu in muList:
	mu, sigma, occupation = findMu(mu_a,mu_b,n,1e-5,Betas[i],U,Ef,p1,nkDMFT,nk,nvDMFT,nv,DMFTiter)

	print(mu, sigma, occupation)
	
	muList[i] = mu
	
muList.tofile("Mu_U05")

	
	
	
	
	
	
	

