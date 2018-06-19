import subprocess
import numpy as np
from shutil import copyfile
import math

'''finding the chemical potential mu for given parameters so that occupation of mobile electrons is n'''
def findMu(mu_a,mu_b,n,error,beta,U,Ef,p1,nkDMFT,nk,nvDMFT,nv,DMFTiter):
	
	#start intervall
	a=mu_a
	b=mu_b
	
	#wanted deviation (sigma) from half filling
	sigmaWanted = n-0.5
	
	#maximal number of steps
	maxSteps = 1000
	
	#initialize Sigma Corrections to zero
	dumsigcor = np.zeros(( nk*nk*2*nv ) , dtype = complex)
	dumsigcor.tofile("SigmaCork")

	#bisection
	for i in range(0,maxSteps):
	
		mu = a + (b - a)/2.



		#call DMFT
		subprocess.call( ["./SCDMFT.out" , str (beta) , str (U) , str (mu) , str (Ef) , str(p1) , str (nkDMFT) , str (nk) , str (nvDMFT) , str (nv) , str(DMFTiter)] )

		G1 = np.fromfile("G1", dtype = complex)

		#calculate deviation from half filling -> sigma
		sigma = sum(G1).real/beta
		occupation = sigma + 0.5
		
		if (sigma - sigmaWanted) >= 0:
			a = mu
		else:
			b = mu
			
		print("T = " + str(1/beta))
		print(mu, sigma)
			
		if ((sigma - sigmaWanted) < error) and ((sigma - sigmaWanted) > -error):
			print("-------------------------")
			print("T = " + str(1/beta) + ": FOUND")
			print("-------------------------")
			return mu, sigma, occupation

		
	
'''calculating the corrected Greens function of DMFT using corrections of DF parquet calculations'''
def calcGCor(pathG1, pathSigma, pathSigmaCork, nk, nv, nvDMFT, beta, mu):
	pi = math.pi
	pistep = 2*pi/float(nk)	

	#read in the DMFT quantities and the corrections
	Gloc = np.fromfile(pathG1, dtype = complex)
	SigmaDMFT = np.fromfile(pathSigma, dtype = complex)
	SigmaCorr = np.fromfile(pathSigmaCork, dtype = complex)
	SigmaCorr = SigmaCorr.reshape((nk,nk,2*nv))

	#initialize quantities
	Epsilon = np.zeros((nk,nk), dtype = float)
	iNu = np.zeros((2*nvDMFT), dtype = complex)
	GlocNew = np.zeros((2*nvDMFT), dtype = complex)
	GDMFT = np.zeros((nk,nk,2*nvDMFT), dtype = complex)
	Gnew = np.zeros((nk,nk,2*nvDMFT), dtype = complex)

	#calculate energy epsilon
	for kx in range(0,nk):
		epsx = math.cos(kx*pistep-pi)

		for ky in range(0,nk):
			Epsilon[kx][ky] = -0.5*(epsx + math.cos(ky*pistep-pi))

	#calculate Matsubara frequencies
	for v in range(0,2*nvDMFT):
		iNu[v] = 1j*(2*(v-nvDMFT)+1)*pi/beta

	#calculate GDMFT and Gnew
	for v in range(0,2*nvDMFT):
		for kx in range(0,nk):
			for ky in range(0,nk):
				GDMFT[kx][ky][v] = 1./(iNu[v]-Epsilon[kx][ky]+mu - SigmaDMFT[v])
				if (v >= (nvDMFT-nv)) and (v < (nvDMFT + nv)):
					dfv = v - (nvDMFT-nv)
					Gnew[kx][ky][v] = 1./(iNu[v]-Epsilon[kx][ky]+mu - SigmaDMFT[v] - SigmaCorr[kx][ky][dfv])
				else:
					Gnew[kx][ky][v] = GDMFT[kx][ky][v]

	#calculate new Gloc
	for v in range(0,2*nvDMFT):
		Gsum = 0
		for kx in range(0,nk):
			for ky in range(0,nk):
				Gsum += Gnew[kx][ky][v]
		GlocNew[v] = Gsum/float(nk*nk)	


	#GDMFT.tofile("GDMFT")
	#Gnew.tofile("Gnew")
	GlocNew.tofile("GlocNew")

	
	
	
	
	
	

