import subprocess
import numpy as np
from shutil import copyfile
import os
'''
Calls SCDMFT.out to calculate the DMFT quantities and, starting with 0 as dual self energy Sigma,
calls DFParquet.out to calculate the DF quantities iteratively. In each step, DualSig and
Gdual are copied to the directories ./SigmaDual and ./GDual.
'''
Ef = 1.
p1 = 0.5
nk = 6
nkDMFT = 128
nv = 20#128
nvDMFT = 100#128
DMFTiter = 40
relconv = 0.5
DFiter = 12
iteration = 3 #iteration steps

#basepath = os.getcwd()

#letters = ('A','B','C','D','E')
#numbers = ('1','2','3','4','5')

#for i in range(5):
  #for j in range(5):
    #os.mkdir(letters[i] + numbers[j])

#os.chdir(basepath + "/Workplace")

copfiles = ("G1","LambdaUp","Sigma")


for i in range(1):
  for j in range(1):
    
    T = 0.02 * (i+1)
    U = 0.2 * (j+1)
    mu = .5 * U
    beta =1./T
    
    #initialize Sigma Corrections to zero
    dumsigcor = np.zeros(( nk*nk*2*nv ) , dtype = complex)
    dumsigcor.tofile("SigmaCork")
    
    #call DMFT
    subprocess.call( ["./SCDMFT.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nkDMFT) , str (nk) , str (nvDMFT) , str (nv) , str(DMFTiter)] )
    
    #for fil in copfiles:
      #copyfile(fil,"../" + letters[i] + numbers[j] + "/" + fil)
    #initialize Dual Sigma to zero
    dumsigdual = np.zeros(( nk*nk*2*nv ) , dtype = complex)
    dumsigdual.tofile("DualSig")
    
    
    for it in range(0,iteration):
      
      #print("============================================")
      #print("Starting DF iteration " + str(i) + "/" + str(iteration-1))
      #print("============================================")
      
      #call DFParquet
      subprocess.call( ["./DFphLadder.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nk) , str (nkDMFT) , str (nv) , str (nvDMFT) , str(relconv) , str(DFiter)] )
      subprocess.call( ["./DFphLadder.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nk) , str (nkDMFT) , str (nv) , str (nvDMFT) , str(relconv) , str(DFiter)] )
      
      #copy Dual Sigma and G Dual
      #directories SigmaDual and GDual have to be already there!




