import subprocess
import numpy as np
import matplotlib.pyplot as plt
from shutil import copyfile

'''
Calculate self energy and Green's function with a Dual Fermion parquet iteration, including
one-particle reducible contributions
'''

beta = 12.5
U = 1.
mu = .5
Ef = 1.
p1 = 0.5
nk = 6
nkDMFT = 128
nv = 20#128
nvDMFT = 20#128
DMFTiter = 40
relconv = 0.5
DFiter = 5

iteration = 5 #iteration steps

#initialize Sigma Corrections to zero
dumsigcor = np.zeros(( nk*nk*2*nv ) , dtype = complex)
dumsigcor.tofile("SigmaCork")

#call DMFT
subprocess.call( ["./SCDMFT.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nkDMFT) , str (nk) , str (nvDMFT) , str (nv) , str(DMFTiter)] )

#initialize Dual Sigma to zero
dumsigdual = np.zeros(( nk*nk*2*nv ) , dtype = complex)
dumsigdual.tofile("DualSig")

for i in range(0,iteration):
	
	print("============================================")
	print("Starting iteration " + str(i) + "/" + str(iteration-1))
	print("============================================")
	
	#call DFParquet
	subprocess.call( ["./DFParquet.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nk) , str (nkDMFT) , str (nv) , str (nvDMFT) , str(relconv) , str(DFiter)] )
	
	#read SigmaDual
	SigmaDual = np.fromfile("DualSig", dtype = complex)
	SigmaDual = SigmaDual.reshape((nk,nk,2*nv))
	
	#read Gloc
	Gloc = np.fromfile("G1", dtype = complex)

	#calculate mapped SigmaDual
	SigmaDualLadder = np.zeros((nk,nk,2*nv), dtype = complex)
	for k in range(0,nk):
		for j in range(0,nk):
			for v in range(0,2*nv):
				SigmaDualLadder[k][j][v]=SigmaDual[k][j][v]/(1-Gloc[v]*SigmaDual[k][j][v])
				
	SigmaDualLadder.tofile("DualSig")
	
	#calculate Sigma corrections
	SigmaCorr = np.zeros((nk,nk,2*nv), dtype = complex)
	for k in range(0,nk):
		for j in range(0,nk):
			for v in range(0,2*nv):
				SigmaCorr[k][j][v]=SigmaDualLadder[k][j][v]/(1+Gloc[v]*SigmaDualLadder[k][j][v])
				
	SigmaCorr.tofile("SigmaCorr")
	
	#copy Dual Sigma, G Dual and Sigma corrections
	copyfile("DualSig","./SigmaDualMap/DualSigMap" + str(i))
	copyfile("Gdual","./GDualMap/GdualMap" + str(i))
	copyfile("SigmaCorr", "./SigmaCorrMap/SigmaCorr" + str(i))
	