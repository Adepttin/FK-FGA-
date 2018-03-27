import matplotlib.pyplot as plt
import numpy as np
import cmath as cm
import subprocess
from SecondOrderFunctions import *

#Paramters
#nk=16
nkList=[4,8,16,32,64]
nv=20
beta=12.5

#DMFT Parameters
nk = 6
U = 1 #1.
mu = 0.5 #.5
Ef = 1. #1.
p1 = 0.5 #0.5
nkDMFT = 128
nvDMFT = 100#128
DMFTiter = 40

#for nk in nkList:
print("---------------------------")
print("Lattice: " + str(nk) + "x" + str(nk))
print("---------------------------")

#initialize Sigma Corrections to zero
dumsigcor = np.zeros(( nk*nk*2*nv ) , dtype = complex)
dumsigcor.tofile("SigmaCork")

#call DMFT
subprocess.call( ["./SCDMFT.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nkDMFT) , str (nk) , str (nvDMFT) , str (nv) , str(DMFTiter)] )

#exponential coefficients
exponent, exponent_neg = setUpExponent(nk)

#reading in Floc and G0dual
print("Reading in Floc and G0dual")

F = readFloc(nv)
#print(len(F))
#plt.imshow(F.real, cmap='hot', interpolation='nearest')
#plt.colorbar()
#plt.show()

Gdual = readGdual(nk,nv)
#plt.imshow(Gdual[:,:,10].real, cmap='hot', interpolation='nearest')
#plt.colorbar()
#plt.show()

#Fourier Transformation of Gdual
print("Fourier Transform G0dual")

#for positive r
Gx = np.zeros((nk,nk,2*nv), dtype = complex)
fourierTransform(Gx, Gdual, nk, nv, exponent)

#for negative r
Gx_neg = np.zeros((nk,nk,2*nv), dtype = complex)
fourierTransform(Gx_neg, Gdual, nk, nv, exponent_neg)

#calculate Sigma in real space
print("Calculate Sigma")
Sigmax = np.zeros((nk,nk,2*nv), dtype = complex)
calculateSigma(Sigmax, F, Gx, Gx_neg, nk, nv, beta)

#Fourier Transformation of Sigma back to momentum space
print("Fourier Transformation back")
Sigmak = np.zeros((nk,nk,2*nv), dtype = complex)
fourierTransform(Sigmak, Sigmax, nk, nv, exponent_neg, back = True)

#write to "DualSigma"
Sigmak.tofile("SecondSigma" + str(nk))

