using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <complex>
#include <algorithm>
#include "routines.cpp"

//class for DF pp ladder calculation based on DMFT results

//--------------------------------------
//routines for calculation of the ladder
//--------------------------------------

//calculate Chiq
int ChiCalc(dcomp ** const Chiq, dcomp * const Gk, const int nk, const int nv, const int* const qtok, const int* const * const kmap, const int* const * const ksum, const int* const * const kdif)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	const int nnk = nk*nk;
	
	const  double norm = nnk;
	int i,j,k,l;
	int k1;
	
	//initialize Chiq to zero
	for(i=0;i<ndistk;i++) //q
	{
		for(k=-nv;k<nv;k++) //v
		{
			for(l=-nv;l<nv;l++) //v'
			{
				Chiq[i][k*2*nv+l] = 0.;
			}
		}
	}
	
	//calculate Chiq
	for(i=0;i<ndistk;i++) //q
	{
		for(j=0;j<nnk;j++) //k
		{
			//q - k
			k1 = kdif[qtok[i]][j];
			
			for(k=-nv;k<nv;k++) //v
			{
				for(l=-nv;l<=k;l++) //v'
				{
					Chiq[i][k*2*nv+l] += (Gk[j*2*nv + k] * Gk[k1*2*nv + l])/norm;
				}
			}
		}
	}
	
	//expand Chiq to full frequency range
	for(i=0;i<ndistk;i++) //q
	{
		for(k=-nv;k<nv;k++) //v
		{
			for(l=k;l<nv;l++) //v'
			{
				Chiq[i][k*2*nv+l] = Chiq[i][l*2*nv+k];
			}
		}
	}
	
	return(0);
}

//calculate vertex function F as given in pp ladder
int FqCalc(dcomp ** const Fq, dcomp ** const Chiq, dcomp * const a, const int nk, const int nv)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	
	int i,k,l;
	
	for(i=0;i<ndistk;i++) //q
	{
		for(k=-nv;k<nv;k++) //v
		{
			for(l=-nv;l<nv;l++) //v'
			{
				Fq[i][k*2*nv+l] = a[k]*a[l] / (1. - a[k]*a[l]*Chiq[i][k*2*nv+l]);
			}
		}
	}
	
	return(0);
}

//--------------------------------------------------------------------
//initialisation of k quantities, including some symmetry maps for BZ
//--------------------------------------------------------------------

//initialise ksum
//ksum gives the sum of two momenta
int initksum(int*const*const ksum , const int nk)
{
	int kx1,kx2,ky1,ky2;
	int sum1,sum2;
	
	for(kx1 = 0; kx1 < nk; kx1++)
	{
		for(kx2 = 0; kx2 < nk; kx2++)
		{
			sum1 = (kx1+kx2+nk/2)%nk;
			for(ky1 = 0; ky1 < nk; ky1++)
			{
				for(ky2 = 0; ky2 < nk; ky2++)
				{
					sum2 = (ky1+ky2+nk/2)%nk;
					
					ksum[kx1*nk+ky1][kx2*nk+ky2] = (sum1*nk+sum2);
				}
			}
		}
	}
	return(0);
}

//initialise kdif
//kdif gives the difference of two momenta
int initkdif(int*const*const kdif , const int nk)
{
	int kx1,kx2,ky1,ky2;
	int dif1,dif2;
	
	for(kx1 = 0; kx1 < nk; kx1++)
	{
		for(kx2 = 0; kx2 < nk; kx2++)
		{
			dif1 = (kx1-kx2+(3*nk)/2)%nk;
			for(ky1 = 0; ky1 < nk; ky1++)
			{
				for(ky2 = 0; ky2 < nk; ky2++)
				{
					dif2 = (ky1-ky2+(3*nk)/2)%nk;
					
					kdif[kx1*nk+ky1][kx2*nk+ky2] = (dif1*nk+dif2);
				}
			}
		}
	}
	return(0);
}

//auxiliary function for use in initktoq
inline int RC( const int k1 , const int k2)
{
	return ((k1*(k1+1))/2 + k2);
}

//initialise ktoq
//ktoq converts a k index on full BZ to a q index on irreducible BZ
int initktoq(int*const ktoq , const int nk)
{
	int kx,ky;
	int xbas,ybas;
	int a,b;
	
	for(kx=0; kx<nk; kx++)
	{
		xbas = min(kx,(nk-kx));
		
		for(ky=0; ky<nk; ky++)
		{
			ybas = min(ky,(nk-ky));
			
			a = max(xbas,ybas);
			b = min(xbas,ybas);
			
			ktoq[kx*nk+ky] = RC(a,b);
		}
	}
	
	return(0);
}

//iniitalise qtok
//qtok converts a q index on irreducible BZ to a k index on full BZ
int initqtok(int*const qtok , const int nk)
{
	int kx,ky;
	int c1;
	
	for(kx=0; kx<=nk/2; kx++)
	{
		c1 = ((kx+1)*kx)/2;
		
		for(ky=0; ky<=kx; ky++)
		{
			qtok[c1+ky] = kx*nk+ky;
		}
	}
	
	return(0);
}

//initialise ksym
//ksym gives symmetry index for irreducible BZ
int initksym(int*const ksym , const int nk)
{
	int kx,ky;
	
	for(kx=0; kx <= nk/2; kx++)
	{
		for(ky=0; ky<= kx; ky++)
		{
			ksym[kx*nk+ky] = 0;
		}
	}
	
	for(kx=1; kx <= nk/2; kx++)
	{
		for(ky=0; ky < kx; ky++)
		{
			ksym[ky*nk+kx] = 1;
		}
	}
	
	for(kx=1; kx < nk/2; kx++)
	{
		for(ky=0; ky<= kx; ky++)
		{
			ksym[(nk-kx)*nk+ky] = 2;
		}
	}
	
	for(kx=2; kx <= nk/2; kx++)
	{
		for(ky=1; ky < kx; ky++)
		{
			ksym[(nk-ky)*nk+kx] = 3;
		}
	}
	
	for(kx=1; kx < nk/2; kx++)
	{
		for(ky=0; ky <= kx; ky++)
		{
			ksym[ky*nk+(nk-kx)] = 4;
		}
	}
	
	for(kx=2; kx <= nk/2; kx++)
	{
		for(ky=1; ky < kx; ky++)
		{
			ksym[kx*nk+(nk-ky)] = 5;
		}
	}
	
	for(kx=1; kx < nk/2; kx++)
	{
		for(ky=1; ky <= kx; ky++)
		{
			ksym[(nk-ky)*nk+(nk-kx)] = 6;
		}
	}
	
	for(kx=2; kx < nk/2; kx++)
	{
		for(ky=1; ky < kx; ky++)
		{
			ksym[(nk-kx)*nk+(nk-ky)] = 7;
		}
	}
	
	return(0);
}

//initialise kmap
//kmap gives index on irreducible BZ
int initkmap(int*const*const kmap , const int nk)
{
	//Initialize Symmetry k map
	int kx,ky;
	int coord;
	int i,j;
	
	for(kx=0; kx < nk; kx++)
	{
		for(ky=0; ky < nk; ky++)
		{
			i = kx;
			j = ky;
			coord = i*nk + j;
			kmap[0][coord] = coord;
			kmap[1][coord] = j*nk + i;
			kmap[2][coord] = ((nk-i)%nk)*nk + j;
			kmap[3][coord] = j*nk + ((nk-i)%nk);
			kmap[4][coord] = ((nk-j)%nk)*nk + i;
			kmap[5][coord] = i*nk + ((nk-j)%nk);
			kmap[6][coord] = ((nk-j)%nk)*nk + ((nk-i)%nk);
			kmap[7][coord] = ((nk-i)%nk)*nk + ((nk-j)%nk);
		}
	}
	
	return(0);
}

//----------------------------------------------------
//calculation of dual self energy - equation of motion
//----------------------------------------------------

//calculate dual self energy
int SigmaCalc(dcomp * DualSig, dcomp ** const Fq, dcomp * const Gk, const int nk , const int nv, const int* const ktoq, const int* const ksym, const int* const qtok, const int* const * const kmap, const int* const * const ksum, const int* const * const kdif)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	const int nnk = nk*nk;
	
	const  double norm = nnk;
	int j,k,l;
	
	int k1,q1,k2,s1;
	
	//initialise dual self energy to zero
	for(j=0;j<nnk;j++) //k
	{
		for(k=-nv;k<nv;k++) //v
		{
			DualSig[2*nv*j + k] = 0.;
		}
	}
	
	for(k1=0;k1<nnk;k1++) //k
	{		
		for(j=0;j<nnk;j++) //q
		{
			//q as q index on irreducible BZ
			q1 = ktoq[j];
			//k - q
			k2 = kdif[k1][j];
			
			for(k=-nv;k<nv;k++) //v
			{
				for(l=-nv;l<nv;l++) //v'
				{
					DualSig[2*nv*k1 + k] -= Fq[q1][2*nv*k + l] * Gk[k2*2*nv + l]/norm;
				}
				DualSig[2*nv*k1 + k] += Fq[q1][2*nv*k + k] * Gk[k2*2*nv + k]/norm;
			}
		
		}
	}
		
	return(0);
}



class DFppParams
{
	public:
	
	//lattice size
	
	//nk - size of DF lattice
	//nkin - size of DMFT lattice
	//nv - number of DF Matsubara frequencies
	//nvin - number of DMFT Matsubara frequencies
	int nk, nkin, nv, nvin;
	//nnk = nk*nk
	int nnk;
	//number of points in irreducible BZ
	int ndistk;
	
	//one-particle quantities
	
	//dual self energy Sigma(k,v)
	dcomp* Sigmadual;
	//resulting DF corrections Sigma_corr(k,v)
	dcomp* Sigmacor;
	//bare dual propagator G_0(k,v)
	dcomp* G0dual;
	//dual propagator G(k,v)
	dcomp* Gdual;
	//local propagator of real fermions G_loc(v)
	dcomp* Gloc;
	
	//k-grid quantities
	
	//sum of two momenta
	int ** ksum;
	//difference of two momenta
	int ** kdif;
	//symmetry map on irreducible BZ
	int ** kmap;
	//symmetry index for use in kmap
	int * ksym;
	//convert k to q index (from full to irreducible BZ)
	int * ktoq;
	//convert q to k index (from irreducible to full BZ)
	int * qtok;
	
	//Vertex quantities
	dcomp* a;
	dcomp ** Fup;	
	dcomp ** Chiq;
	
	//constructor
	DFppParams( const int innk , const int innv , const int innkin , const int innvin )
	{
		nk = innk;
		nv = innv;
		nkin = innkin;
		nvin = innvin;
		
		nnk = nk*nk;
		ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
		
	}
	
	int InitialiseOneParticle()
	{
		G0dual = new dcomp [2*nv*nnk]; 
		G0dual += nv;
		
		Gdual = new dcomp [2*nv*nnk]; 
		Gdual += nv;
		
		Sigmadual = new dcomp [2*nv*nnk]; 
		Sigmadual += nv;
		
		Sigmacor = new dcomp [2*nv*nnk]; 
		Sigmacor += nv;
		
		Gloc = new dcomp [2*nv]; 
		Gloc += nv;
		
		return(0);
	}
	
	int DeleteOneParticle()
	{
		G0dual -= nv;
		delete G0dual;
		
		Gdual -= nv;
		delete Gdual;
		
		Sigmadual -= nv;
		delete Sigmadual;
		
		Sigmacor -= nv;
		delete Sigmacor;
		
		Gloc -= nv;
		delete Gloc;
		
		return(0);
	}
	
	int InitialiseKQuantities()
	{
		
		ksum = new int *[nnk];
		kdif = new int *[nnk];
		kmap = new int *[8];
		ksym = new int [nnk];
		ktoq = new int [nnk];
		qtok = new int [ndistk];
		
		{
			int i;
			for (i = 0; i< 8; i++)
			{
				kmap[i] = new int [nnk];
			}
			for (i = 0; i< nnk; i++)
			{
				ksum[i] = new int [nnk];
				kdif[i] = new int [nnk];
			}
		}
		
		initksum(ksum , nk);
		initkdif(kdif , nk);
		initktoq(ktoq , nk);
		initksym(ksym , nk);
		initkmap(kmap , nk);
		initqtok(qtok , nk);
		
		return(0);
	};
	
	int DeleteKQuantities()
	{
		int i;
		for (i = 0; i< 8; i++)
		{
			delete kmap[i];
		}
		for (i = 0; i< nnk; i++)
		{
			delete ksum[i];
			delete kdif[i];
		}
		
		delete ksum;
		delete kdif;
		delete kmap;
		delete ktoq;
		delete qtok;
		delete ksym;
		
		return(0);
	};
	
	int InitialiseVertexStorage()
	{
		
		int i, j;
		
		a = new dcomp [2*nv];
		a += nv;
		
		Fup = new dcomp *[ndistk];
		
		Chiq = new dcomp *[ndistk];
		
        for (i = 0; i< ndistk; i++)
		{
			*(Fup + i) = new dcomp [4*nv*nv];
			*(Chiq + i) = new dcomp [4*nv*nv];
			
			*(Fup + i) += (2*nv + 1)*nv;
			*(Chiq + i) += (2*nv + 1)*nv;
		}
		
		return(0);
	};
	
	int DeleteVertexStorage()
	{
		int i;
		
		for (i = 0; i< ndistk; i++)
		{
			Fup[i] -= (2*nv + 1)*nv;
			Chiq[i] -= (2*nv + 1)*nv;
			
			delete Fup[i];
			delete Chiq[i];
		}
		
		delete Fup;
		
		delete Chiq;
		
		return(0);
	};
	
	//read in DMFT results
	int ReadDMFT()
	{
		int i, j;
		
		//read in a
		
		dcomp* ain = new dcomp [2*nvin];
		
		readbin ("a" , ain , 2*nvin );
		
		ain += nvin;
		
		for(i=-nv; i<nv;  i++)
		{
			a[i] = ain[i];
		}
		
		{
			int i,j;
			for (i = 0; i< nnk; i++)
			{
              for (j = 0; j< 1; j++)
              {
                  cout << i << "/" << nnk << "    " << j << "/" << nnk << "    " << kdif[i][j] << "/" << nnk << endl;
              }
			}
		}
		cout << "/areset                      /" << endl;
		
		ain -= nvin;
		
		delete ain;
		
		//read in local Green's function
		
		dcomp* Gin = new dcomp [2*nvin];
		
		readbin ("G1" , Gin , 2*nvin );
		
		Gin += nvin;
		for(i=-nv; i<nv;  i++)
		{
			Gloc[i] = Gin[i];
		}
		Gin -= nvin;
		
		delete Gin;
		
		//read in bare dual Green's function
		
		readbin ("G0dual" , G0dual-nv , 2*nv*nnk );
		
		for(i=0; i<nnk;  i++)
		{
			for(j=-nv; j<nv; j++)
			{
				Gdual[i*2*nv + j] = G0dual[i*2*nv + j];
			}
		}
		
		return(0);
	}
	
	//set self energy corrections equal to dual self energy
	int DualToRealSig()
	{
		int i,j;
		for(i=0; i< nnk; i++)
		{
			for(j=-nv; j< nv; j++)
			{
				Sigmacor[i*2*nv + j] = Sigmadual[i*2*nv + j];
			}
		}
		
		return(0);
	}
	
	//use mapping to calculate self energy corrections out of dual self energy
	int MappedDualToRealSig()
	{
		int i,j;
		for(i=0; i< nnk; i++)
		{
			for(j=-nv; j< nv; j++)
			{
				Sigmacor[i*2*nv + j] = Sigmadual[i*2*nv + j]/ (1. +  Sigmadual[i*2*nv + j] * Gloc[j]);
			}
		}
		
		return(0);
	}
	
	//choose which way self energy corrections are calculated
	int FlexDualToRealSig(int mode)
	{
		switch(mode)
		{
			case 0:
				this->DualToRealSig();
				break;
			case 1: 
				this->MappedDualToRealSig();
				break;
			default:
				return (1);
		}
		
		return (0);
	}
	
	//set dual self energy equal to self energy corrections
	int RealToDualSig()
	{
		int i,j;
		for(i=0; i< nnk; i++)
		{
			for(j=-nv; j< nv; j++)
			{
				Sigmadual[i*2*nv + j] = Sigmacor[i*2*nv + j];
			}
		}
		
		return(0);
	}
	
	//use mapping to calculate dual self energy out of self energy corrections
	int MappedRealToDualSig()
	{
		int i,j;
		for(i=0; i< nnk; i++)
		{
			for(j=-nv; j< nv; j++)
			{
				Sigmadual[i*2*nv + j] = Sigmacor[i*2*nv + j]/ (1. -  Sigmacor[i*2*nv + j] * Gloc[j]);
			}
		}
		
		return(0);
	}
	
	//choose which way dual self energy is calculated
	int FlexRealToDualSig(int mode)
	{
		switch(mode)
		{
			case 0:
				this->RealToDualSig();
				break;
			case 1: 
				this->MappedRealToDualSig();
				break;
			default:
				return (1);
		}
		
		return (0);
	}
	
	//calculate Green's function out of dual self energy via Dyson equation
	int UpdateGdual()
	{
		int i,j;
		for(i=0; i< nnk; i++)
		{
			for(j=-nv; j< nv; j++)
			{
				Gdual[i*2*nv + j] = G0dual[i*2*nv + j]/ (1. -  Sigmadual[i*2*nv + j] * G0dual[i*2*nv + j]);
			}
		}
		
		return(0);
	}
	
	//read in self energy corrections
	int ReadSigCors()
	{
		readbin ("SigmaCork" , Sigmacor-nv , 2*nv*nk*nk );
		
		return(0);
	}
	
	//write self energy corrections
	void WriteSigCors()
	{
		writebin ("SigmaCork" , Sigmacor-nv , 2*nv*nnk );
	}
	
	//read in dual self energy
	int ReadDualSig()
	{
		readbin("DualSig", Sigmadual-nv, 2*nv*nnk);
		
		return(0);
	}
	
	//write dual self energy
	void WriteDualSig()
	{
		writebin ("DualSig" , Sigmadual-nv , 2*nv*nnk );
	}
	
	//write dual Green's function
	void WriteGdual()
	{
		writebin ("Gdual" , Gdual-nv , 2*nv*nnk );
	}
	
	//call pp ladder calculation
	int LadderCalc()
	{
		ChiCalc(Chiq, Gdual, nk, nv, qtok, kmap, ksum, kdif);
		FqCalc(Fup, Chiq, a, nk, nv);
		
		return(0);
	}
	
	//calculate dual self energy
	int SigCalc()
	{
		SigmaCalc(Sigmadual, Fup, Gdual, nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);
		
		return(0);
	}
};

//class for calculating the current-current correlation function

//prepare velocity-x-operator for nk*nk k-mesh
//nearest neighbour-hopping square lattice is hardcoded so far
template <typename numbertype>
int calcvx(const int nk, numbertype* const vx)
{
	int i,j;
	numbertype v;
	numbertype pi = acos( (-1.));
	numbertype pistep = 2.*pi / nk;
	const double t = 0.25;
	
	for(i = 0; i < nk; i++)
	{
		v = 2. * t * sin( i * pistep - pi );
		for(j = 0; j < nk; j++)
		{
			*(vx + i*nk + j) = v;
		}
	}
	
	return(0);
}

//prepare fermionic Matsubara frequencies
template <typename cpltype, typename fltype>
int PrecalciNu( cpltype* const iNu , const fltype beta, int nv)
//precalculates a vector of Matsubara frequencies for use in some functions.
{
	int i;
	fltype pibeta = acos( fltype(-1.))/beta;
	
	for(i=-nv;i< nv;i++)
	{
		*(iNu + i) = cpltype( 0. , pibeta * (2*i + 1) );
	}
	
	return (0);
}

//calculate bubble contribution
//use Gk (already dressed by non-local self energy corrections), Matsubara frequencies and x-velocities vx as input
//correct for finite frequency box by subtracting the trivial propagators 1/(iv * (iv + iw)) which give 0.25 for w = 0 and 0 otherwise when summed over frequencies
//write result to ohm
template <typename numbertype , typename fltype>
int calcohm(const numbertype* const Gk, const numbertype* const Matsus, const fltype* const vx, const int nv2, const int nw2, const int nk, const fltype beta, numbertype* const ohm )
{
	int i,j,k,l;
	numbertype dummy, dumg;
	const int norm = nk*nk;
	numbertype duma, dumb;
	const fltype flnorm = norm;
	
	fltype fldummy, zeroohm;
	
	//auxiliary term
	
	zeroohm = 0.;
	
	for(j = 0; j < nk; j++)
	{
		fldummy = 0.;
		
		for(l = 0; l < nk; l++)
		{
			fldummy += (-0.25) * *(vx + j*nk + l) * *(vx + j*nk + l);
			
		}
		zeroohm += fldummy;
	}
	
	zeroohm = zeroohm/flnorm;
	
	//calculate bubble	
	
	for(k = 0; k < nw2+1; k++) //w
	{
		//initialise bubble to zero
		*(ohm + k) = 0.;
		
		//dummy is a primitive anti-absorbtion measure for numerical summation
		for(j = 0; j < nk; j++) //k_x
		{
			dummy = 0.;
			
			for(l = 0; l < nk; l++) //k_y
			{
				for(i = -nv2; i < nv2-k; i++)
				{
					//G(k,v)
					duma = *(Gk + (j*nk + l)*2*nv2 + i);
					//G(k,v+w)
					dumb = *(Gk + (j*nk + l)*2*nv2 + i+k);
					//dummy propagator
					dumg = 1./(Matsus[i] * Matsus[i+k]);
					
					//subtract dummy propagator in sum
					dummy += ( duma * dumb - dumg ) * *(vx + j*nk + l) * *(vx + j*nk + l);
				}
				
			}
			
			*(ohm + k) += dummy;
		}
		
		//normalization
		*(ohm + k) = *(ohm + k) / (flnorm);
		*(ohm + k) = *(ohm + k) / (beta*beta);
	}
	
	//add auxiliary term to account for subtraction of dummy propagator in sum
	*(ohm) += zeroohm;
	
	//expand to negative Matsubara frequencies
	for(k = 1; k < nw2+1; k++)
	{
		*(ohm - k) = *(ohm + k);
	}
		
	return(0);
}


//calculate vertex corrections
//use Fq, Gk (already dressed by non-local self energy corrections), Matsubara frequencies and x-velocities vx as input
//correct for finite frequency box by subtracting the trivial propagators C * 1/(iv * (iv + iw))^2 which give w-constant contributions per kk' set
//write result to ohm
template <typename numbertype , typename fltype>
int calcconohm(const numbertype* const* const Fqdown, const numbertype* const Gk, const numbertype* const Matsus, const fltype* const vx, const int nv2, const int nw2, const int nk, const fltype beta, numbertype* const ohm , const int *const ktoq, const int* const * const kdif)
{
	int i,k,kp,qind,qredind,v;
	numbertype dummya, dummyb, dummyc, dumg, Casym;
	numbertype dumg2;
	fltype fldummy;
	numbertype * wref;
	const int norm = nk*nk;
	const fltype flnorm = norm;
	const fltype pi = acos( (-1.));
	
	wref = new numbertype[nw2 + 1];
	
	*wref = beta*beta / 48.;
	
	fldummy = beta*beta/(8. * pi * pi);
	for(i = 1; i < nw2 + 1; i++)
	{
		*(wref+i) = fldummy / (i*i);
	}
	
	for(i = 0; i < nw2 + 1; i++) // w
	{
		*(ohm + i) = 0.;
		
		for(qind = 0; qind < norm; qind++) //q
		{
			//q as q index on irreducible BZ
			qredind = ktoq[qind];
			
			dummya = 0.;
			
			for(k = 0; k < norm; k++) //k
			{
				//q - k
				kp = kdif[qind][k];
				
				dummyb = 0.;
				
				Casym = Fqdown[qredind][(nv2-1)*(2*nv2)-nv2];
				
				for(v = -nv2; v < nv2-i; v++) //v
				{
					//dummy propagator for substraction
					dumg = 1./(Matsus[v] * Matsus[i+v]); 
					dumg = dumg*dumg;
					//G(k,v)*G(q-k,v)*G(q-k,v + w)*G(k,v + w)
					dumg2 = Gk[k*2*nv2 + v] * Gk[kp*2*nv2 + v] * Gk[kp*2*nv2 + (v+i)] * Gk[k*2*nv2 + (v+i)];
					
					//calculate vertex corrections with substraction of dummy term
					dummyb += ( (Fqdown[qredind][(v)*(2*nv2+1)+i] * dumg2) - (Casym * dumg) ) * vx[k] * vx[kp];
				}
				
				//add dummy term again
				dummya += (dummyb/(beta*beta) + Casym * wref[i] * vx[k] * vx[kp]);
			}
			
			//different sign and normalization
			*(ohm + i) -= dummya/(flnorm*flnorm);
		}
	}
	
	//extend to negative Matsubara frequencies
	for(k = 1; k < nw2+1; k++)
	{
		*(ohm - k) = *(ohm + k);
	}
	
	delete wref;
	return(0);
}

class ConductivityObject
{
	public:
	
	//parameters
	
	//nk - lattice size
	//nv - number of Matsubara frequencies
	//nvin - number of DMFT Matsubara frequencies
	int nk, nv, nvin;
	double beta,mu;
	//velocity operator in x-direction
	double* vx;
	//dispersion relation
	double* Ek;
	
	//one-particle quantities
	
	//local self energy
	dcomp* Sigmaloc;
	//self energy corrections
	dcomp* Sigmacor;
	//Matsubara frequencies
	dcomp* Matsus;
	//lattice Green's function
	dcomp* Gk;
	
	//k-quantities
    int ** kdif;
    int * ktoq;
    
	//bubble current-current correlation function
	dcomp* ohmbubble;
	//vertex corrections to current-current correlation function
	dcomp* ohmvertex;
	
	//vertex quantities
	dcomp ** Fdown;
	
	//constructor using ladder class
	template <class Parquet>
	ConductivityObject(Parquet MyParquet, double inbeta, double inmu)
	{
		nk = MyParquet.nk;
		nv = MyParquet.nv;
		nvin = MyParquet.nvin;
		
		Sigmacor = MyParquet.Sigmacor;
		Fdown = (MyParquet.Fup);
		
        kdif = MyParquet.kdif;
        ktoq = MyParquet.ktoq;
        
		beta = inbeta;
		mu = inmu;
	}
	
	int InitialiseStorage()
	{
		vx = new double[nk*nk];
		Ek = new double[nk*nk];
		
		Sigmaloc = new dcomp[2*nvin];
		Sigmaloc = Sigmaloc + nvin;
		
		Gk = new dcomp[2*nv*nk*nk];
		Gk = Gk + nv;
		
		Matsus = new dcomp[2*nv];
		Matsus = Matsus + nv;
		
		ohmbubble = new dcomp[2*nv + 1];
		ohmbubble = ohmbubble + nv;
		
		ohmvertex = new dcomp[2*nv + 1];
		ohmvertex = ohmvertex + nv;
		
		return(0);
	}
	
	int DeleteStorage()
	{
		delete vx;
		delete Ek;
		
		Sigmaloc = Sigmaloc - nvin;
		delete Sigmaloc;
		
		Gk = Gk - nv;
		delete Gk;
		
		Matsus = Matsus - nv;
		delete Matsus;
		
		ohmbubble = ohmbubble - nv;
		delete ohmbubble;
		
		ohmvertex = ohmvertex - nv;
		delete ohmvertex;
		
		return(0);
	}
	
	int InitialiseQuantities()
	{
		calcvx(nk, vx);
		calcEk(nk, Ek);
		PrecalciNu( Matsus , beta, nv);
		
		readbin ("Sigma", Sigmaloc-nvin, 2*nvin);
		
		calcGreal( Sigmacor, Sigmaloc, Ek, nv, nk, mu, beta, Gk);

		return(0);
	}
	
	int CalcCondBubble()
	{
		calcohm(Gk, Matsus, vx, nv, nv, nk, beta, ohmbubble);
		return(0);
	}
	
	int CalcCondVertex()
	{
		calcconohm(Fdown, Gk, Matsus, vx, nv, nv, nk, beta, ohmvertex, ktoq, kdif);
		return(0);
	}
	
	int WriteConductivities()
	{
		writebin("BubbleOhm", ohmbubble-nv, 2*nv+1);
		writebin("VertexOhm", ohmvertex-nv, 2*nv+1);
		
		return(0);
	}
	
};

/* not used in the code so far

int SetupQK(int ** const QK , const int nk)
{
 	const int nnk = nk*nk;
 	const int ndistk = (nk/2 + 1)*(nk/2+2)/2;
	
	int i,j,k,l;
	int ca, cb;
	int dumi, dumj;
	
	for(i=0; i<=nk/2;  i++)
	{
		for(j=0; j<=i;  j++)
		{
			ca = RC (i,j);
			for(k=0; k < nk;  k++)
			{
				dumi = ((i+k+nk/2)%nk)*nk;
				cb = k * nk;
				for(l=0; l < nk;  l++)
				{
					dumj = ((j+l+nk/2)%nk);
					
					QK[ca][cb+l] = dumi + dumj;
				}
			}
		}
	}
	
	
	return (0);
} */
