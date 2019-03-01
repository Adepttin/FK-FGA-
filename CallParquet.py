import subprocess
import numpy as np
from shutil import copyfile

#parameters

#inverse temperature
beta = 10. 
#Coulomb interaction strength
U = 0.5
#chemical potential
mu = 0.25
#occupation of localized f-electrons
p1 = 0.5
#lattice size for DF
nk = 6
#lattice size for DMFT
nkDMFT = 128
#number of Matsubara frequencies for DF
nv = 20
#number of Matsubara frequencies for DMFT
nvDMFT = 100
#number of iterations in the self consistent DMFT cycle
DMFTiter = 40
n#umber of Bethe-Salpeter and parquet iterations in DFParquet
DFiter = 7
#number of iterations for the whole DF parquet calculation
iteration = 5

#initialize self energy corrections to zero
dumsigcor = np.zeros((nk*nk*2*nv) , dtype = complex)
dumsigcor.tofile("SigmaCork")

#call SCDMFT
subprocess.call( ["./SCDMFT.out" , str (beta) , str (U) , str (mu) , str(p1) , str (nkDMFT) , str (nk) , str (nvDMFT) , str (nv) , str(DMFTiter)] )

#initialize dual self energy to zero
dumsigdual = np.zeros((nk*nk*2*nv) , dtype = complex)
dumsigdual.tofile("DualSig")

for i in range(0,iteration):

	print("============================================")
	print("Starting DF iteration " + str(i+1) + "/" + str(iteration))
	print("beta = " + str(beta))
	print("============================================")

	#call DFParquet
	subprocess.call( ["./DFParquet.out" , str (beta) , str (U) , str (mu) , str(p1) , str (nk) , str (nkDMFT) , str (nv) , str (nvDMFT) , str(DFiter)] )

	#copy quantities
	copyfile("DualSig","./SigmaDual/DualSig" + str(i))
	copyfile("Gdual","./GDual/Gdual" + str(i))
	copyfile("BubbleOhm","./Ohm/BubbleOhm" + str(i))
	copyfile("VertexOhm","./Ohm/VertexOhm" + str(i))
	copyfile("ChiBubble","./Chi/ChiBubble" + str(i))
	copyfile("ChiVertex","./Chi/ChiVertex" + str(i))
