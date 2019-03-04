'''
DF ladder calculation for the FKM.
Calls SCDMFT.out to calculate the DMFT quantities and calls either
DFphLadder.out or DFppLadder.out to do a ph or pp ladder calculation. 
'''

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
#lattice size for DF ladder
nk = 32
#lattice size for DMFT
nkDMFT = 128
#number of Matsubara frequencies for DF ladder
nv = 40
#number of Matsubara frequencies for DMFT
nvDMFT = 100
#number of iterations in the self consistent DMFT cycle
DMFTiter = 40
#number of iterations for a self consistent DF ladder calculation
iteration = 5

#ph or pp ladder
ladder = "ph"
    
#initialize self energy corrections to zero
dumsigcor = np.zeros(( nk*nk*2*nv ) , dtype = complex)
dumsigcor.tofile("SigmaCork")
    
#call SCDMFT
subprocess.call( ["./SCDMFT.out" , str(beta), str(U), str(mu), str(p1), str(nkDMFT), str(nk), str(nvDMFT), str(nv), str(DMFTiter)] )

#initialize dual self energy to zero
dumsigdual = np.zeros(( nk*nk*2*nv ) , dtype = complex)
dumsigdual.tofile("DualSig")
    
    
for i in range(0,iteration):

    print("============================================")
    print("Starting DF ladder iteration " + str(i) + "/" + str(iteration-1))
    print("beta = " + str(beta))
    print("============================================")

    #call DFphLadder or DFppLadder
    if ladder == "ph":
        subprocess.call( ["./DFphLadder.out" , str(beta), str(U), str(mu), str(p1), str(nk), str(nkDMFT), str(nv), str(nvDMFT)] )
    elif ladder == "pp":
        subprocess.call( ["./DFppLadder.out" , str(beta), str(U), str(mu), str(p1), str(nk), str(nkDMFT), str(nv), str(nvDMFT)] )

    #copy quantities
    copyfile("DualSig","./SigmaDual/DualSig" + str(i))
    copyfile("Gdual","./GDual/Gdual" + str(i))
    copyfile("BubbleOhm","./Ohm/BubbleOhm" + str(i))
    copyfile("VertexOhm","./Ohm/VertexOhm" + str(i))




