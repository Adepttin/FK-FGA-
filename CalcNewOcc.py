import subprocess
import numpy as np
from shutil import copyfile
from routines import *

'''
Calculate new occupation of mobile c electrons resulting from DF corrections to DMFT
'''

#inverse temperature
beta = 10.
#chemical potential
mu = 0.25
#DF lattice size
nk = 6
#number of DF Matsubara frequencies
nv = 20
#number of DMFT Matsubara frequencies
nvDMFT = 100

#path of local DMFT Green's function and local DMFT self energy
pathG1 = "G1"
pathSigma = "Sigma"

#calculate occupation in DMFT
G1 = np.fromfile(pathG1, dtype=complex)
occ_DMFT = sum(G1).real/beta + 0.5

#calculate occupation for DF parquet

#path of self energy corrections (dual self energy)
pathSigmaCork = "DualSig"

#calculate corrected DF Green's function
calcGCor(pathG1, pathSigma, pathSigmaCork, nk, nv, nvDMFT, beta, mu)

#read in resulting corrected local Green's function (DMFT + DF)
GlocNew = np.fromfile("GlocNew", dtype = complex)

#calculate new DF occupation
occupation = sum(GlocNew).real/beta + 0.5

#print new occupation
print(occupation)
	
	
	

	
	
	

