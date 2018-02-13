import subprocess
import numpy as np

beta = 20.
U = 1.
mu = .5
Ef = 1.
p1 = 0.5
nk = 8
nkDMFT = 128
nv = 10#128
nvDMFT = 100#128
DMFTiter = 40

dumsigcor = np.zeros(( nk*nk*2*nv ) , dtype = complex)
dumsigcor.tofile("SigmaCork")

subprocess.call( ["./SCDMFT.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nkDMFT) , str (nk) , str (nvDMFT) , str (nv) , str(DMFTiter)] )
