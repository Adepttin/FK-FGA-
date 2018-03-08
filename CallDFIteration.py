import subprocess
import numpy as np
from shutil import copyfile

beta = 12.5
U = 1.
mu = .5
Ef = 1.
p1 = 0.5
nk = 6
nkDMFT = 128
nv = 20#128
nvDMFT = 100#128
DMFTiter = 40
relconv = 0.5
DFiter = 5

iteration = 10 #iteration steps

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
	print("Starting DF iteration " + str(i) + "/" + str(iteration-1))
	print("============================================")
	
	#call DFParquet
	subprocess.call( ["./DFParquet.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nk) , str (nkDMFT) , str (nv) , str (nvDMFT) , str(relconv) , str(DFiter)] )
	
	#copy Dual Sigma and G Dual
	#directories SigmaDual and GDual have to be already there!
	copyfile("DualSig","./SigmaDual/DualSig" + str(i))
	copyfile("Gdual","./GDual/Gdual" + str(i))
	