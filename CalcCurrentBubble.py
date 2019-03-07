import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import BubbleCurrent as bc

'''
Calculate the bubble term of the current-current correlator based on a 
Green's function including DF corrections by Fourier transformation to
imaginary times
'''

#inverse temperature
beta = 10.
#chemical potential
mu = 0.25
#lattice size
nk = 6
#number of Matsubara frequencies for the bubble calculation
nv = 100
#number of Matsubara frequencies in the DF calculation
nvdf = 20


#read in DF corrected Green's function Gdf
Gdf = np.fromfile("Gnew", dtype=complex)
Gdf = Gdf.reshape((nk,nk,2*nv))

#create bubble object
BubbleObj = bc.BubbleCurrent(nk,nv,beta,mu)

#calculate chi_jj bubble
Chiw = BubbleObj.calcBubble(Gdf)
#take only DF frequency range
Chiw = Chiw[(nv-nvdf):(nv+nvdf+1)]

#save result
Chiw.tofile("OhmBubble")
