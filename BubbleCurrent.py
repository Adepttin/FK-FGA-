import numpy as np

class BubbleCurrent:
	"""class to calculate the bubble term of the current-current correlator"""
		
	def __init__(self,nk,nv,beta,mu):
		"""constructor"""
		
		#parameters
		self.nk = nk #lattice size
		self.nv = nv #number of Matsubara frequencies
		self.beta = beta #inverse temperature
		self.mu = mu #chemical potential
		
		#Matsubara frequencies
		xF = np.arange(-nv,nv,1)
		xB = np.arange(-nv,nv+1,1)
		self.MatsusF = (2*xF+1)*np.pi/beta #fermionic
		self.MatsusB = (2*xB)*np.pi/beta #bosonic
		
		#imaginary times
		xTau = np.arange(0,nv+1,1)
		self.Taus = beta*xTau/float(nv)
		
		#dispersion relation epsilon(k)
		self.tt = 0.25 #hopping parameter
		ks = np.arange(-np.pi, np.pi, 2*np.pi/nk)
		kx,ky = np.meshgrid(ks,ks)
		self.Epsilon = -2*self.tt*(np.cos(kx) + np.cos(ky))
		
		
	def calcBubble(self,Gkv):
		"""calculate Chi_jj bubble"""
		
		#Fourier transform Green's function to imaginary times
		#G(k,nu) --> G(k,tau)
		print("FoTo of Gkv to Gkt")
		Gkt = self.FoToGkv(Gkv)
		
		#calculate bubble term in imaginary times
		print("Calc Chit")
		Chit = self.calcChit(Gkt)
		
		#Fourier transform bubble term back to Matsubara frequencies
		#Chi(tau) --> Chi(omega)
		print("FoTo of Chit to Chiw") 
		Chiw = self.FoToChit(Chit)
		
		return Chiw
		
	def FoToGkv(self,Gkv):
		"""Fourier transform Gkv to Gkt"""
		
		nk = self.nk
		nv = self.nv
		beta = self.beta
		mu = self.mu
		
		#calculate free Green's function
		G0 = self.calcFreeG()
		
		#difference of DF G and free G
		DiffG = Gkv - G0
		
		Gkt = np.zeros((nk,nk,nv+1),dtype=complex)
		
		for kx in range(nk):
			for ky in range(nk):
				for t in range(nv+1):
					fsum = 0
					#Fourier transformation of difference
					for v in range(2*nv):
						fsum += np.exp(-1j*self.MatsusF[v]*self.Taus[t])*DiffG[kx][ky][v]/beta
					#add Fourier transformed G_0 again
					Gkt[kx][ky][t] = fsum - np.exp((beta-self.Taus[t])*(self.Epsilon[kx][ky]-mu))/(1.+np.exp(beta*(self.Epsilon[kx][ky]-mu)))
				
		return Gkt
	
	def calcChit(self,Gkt):
		"""calculate Chi_jj bubble in imaginary times"""
		
		nk = self.nk
		nv = self.nv
		
		Chit = np.zeros((nv+1),dtype=complex)

		for t in range(nv+1):
			tsum = 0
			for kx in range(nk):
				vx = 2*self.tt*np.sin(kx)
				for ky in range(nk):
					tsum += Gkt[kx][ky][t]*Gkt[kx][ky][(nv)-t]*vx*vx
			Chit[t]=tsum/(nk*nk)
			
		return Chit
	
	def FoToChit(self,Chit):
		"""Fourier transform Chi_jj to Matsubary frequencies"""
		
		nv = self.nv
		
		Chiw = np.zeros((2*nv+1),dtype=complex)

		for w in range(2*nv+1):
			wfunc = np.exp(1j*self.MatsusB[w]*self.Taus)*Chit
			Chiw[w] = np.trapz(wfunc,self.Taus)
			
		return Chiw
		
	def calcFreeG(self):
		"""calculate free Green's function"""
		
		nk = self.nk
		nv = self.nv
		
		FreeG = np.zeros((nk,nk,2*nv),dtype=complex)
		for kx in range(nk):
			for ky in range(nk):
				eps = self.Epsilon[kx][ky]
				for v in range(2*nv):
					FreeG[kx][ky][v]=1./(1j*self.MatsusF[v]+self.mu-eps)
					
		return FreeG
	

		
