import subprocess
import numpy as np
from shutil import copyfile
from routines import *

'''
Calculate new occupation of mobile electrons resulting from DF corrections to DMFT (Full Parquet, LadderPH or LadderPP)
'''

#beta = 15.
U = 0.5
#mu = -0.3
Ef = 1.
p1 = 0.5
nk = 6
nkDMFT = 128
nv = 20#128
nvDMFT = 200#128
DMFTiter = 40
relconv = 0.5
DFiter = 5

iteration = 5 #iteration steps

#temperatures for which to calculate
temperatures = [0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04]
betas = [1/t for t in temperatures]

#for file handling
betaNames = [10,11,12,14,16,20,25]

#corresponding chemical potentials
mus = [-0.42327881, -0.42108154, -0.41952515, -0.41874695, -0.41888428, -0.42021179, -0.42337036]

typeList = ["Parquet", "LadderPH", "LadderPP"]

occs = open("Occupations.txt", "w")

occs.write("T\t\t\tmu\t\t\t\tn_DMFT\t\tn_Parquet\tn_ph\t\tn_pp\n")

for i in range(0,len(temperatures)):
	
	beta = betas[i]
	mu = mus[i]
	
	#read in DMFT G and Sigma for current beta
	pathG1 = "./U05/DMFT/beta" + str(betaNames[i]) + "/G1"
	pathSigma = "./U05/DMFT/beta" + str(betaNames[i]) + "/Sigma"
	
	#calculate occupation for DMFT
	G1 = np.fromfile(pathG1, dtype=complex)
	occ_DMFT = sum(G1).real/betas[i] + 0.5
		   
   	occList = []
	
	#calculate occupation for DF Parquet, Ladder ph and Ladder pp
	for modus in typeList:
		
		#read in Sigma Corrections (Dual Sigma)
		pathSigmaCork = "./U05/" + modus + "/beta" + str(betaNames[i]) + "/DF4/SigmaDual/DualSig" + str(iteration-1)
		
		calcGCor(pathG1, pathSigma, pathSigmaCork, nk, nv, nvDMFT, beta, mu)
		
		#copyfile("GlocNew", "./U05/" + modus + "/beta" + str(betaNames[i]) + "/DF4/GlocNew")
		
		#read in resulting corrected local G (DMFT + DF)
		GlocNew = np.fromfile("GlocNew", dtype = complex)

		#calculate new occupation
		occupation = sum(GlocNew).real/betas[i] + 0.5
		
		occList.append(occupation)
	
	occs.write("%1.2f\t\t%1.5f\t\t%1.5f\t\t%1.5f\t\t%1.5f\t\t%1.5f\n" % (temperatures[i], mus[i], occ_DMFT, occList[0], occList[1], occList[2]))	
			
occs.close()
	
	
	

	
	
	

