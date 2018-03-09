import matplotlib.pyplot as plt
import numpy as np
import cmath as cm
import subprocess
from SecondOrderFunctions import *

"""
Calculates the second order contribution of the dual self energy.
Calls SDMFT.out.
"""

#Parameter
beta = 12.5 #20
U = 1 #1.
mu = 0.5 #.5
Ef = 1. #1.
p1 = 0.5 #0.5
nk = 16
nkDMFT = 128
nv = 20#128
nvDMFT = 100#128
DMFTiter = 40

#initialize Sigma Corrections to zero
dumsigcor = np.zeros(( nk*nk*2*nv ) , dtype = complex)
dumsigcor.tofile("SigmaCork")

#call DMFT
print("Call DMFT")
subprocess.call( ["./SCDMFT.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nkDMFT) , str (nk) , str (nvDMFT) , str (nv) , str(DMFTiter)] )

#set up exponential coefficients
exponent, exponent_neg = setUpExponent(nk)

#read in Floc and G0dual
print("Read in Floc and G0dual")
F = readFloc(nv)
Gdual = readGdual(nk,nv)

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
fourierTransform(Sigmak, Sigmax, nk, nv, exponent_neg)

#write to "SecondSigma"
Sigmak.tofile("SecondSigma")

