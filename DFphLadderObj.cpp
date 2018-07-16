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

// Self contained pp-ladder DF based on DMFT results


inline int RC( const int k1 , const int k2)
//
{
	return ((k1*(k1+1))/2 + k2);
}

int SetupQK(int ** const QK , const int nk)
{
// 	const int nnk = nk*nk;
// 	const int ndistk = (nk/2 + 1)*(nk/2+2)/2;
	
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
}


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


int initksym(int*const ksym , const int nk)
{
	//Initialize Symmetryindex for irreducible BZ
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

int ChiCalc(dcomp ** const Chiq, dcomp * const Gk, const int nk, const int nv, const int* const qtok, const int* const * const kmap, const int* const * const ksum, const int* const * const kdif)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	const int nnk = nk*nk;
	
	const  double norm = nnk;
	int i,j,k,l;
	int k1;
	
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
	
	
	for(i=0;i<ndistk;i++) //q
	{
		for(j=0;j<nnk;j++) //k
		{
			k1 = ksum[qtok[i]][j];
			
			for(k=-nv;k<nv;k++) //v
			{
				
				for(l=-nv;l<=k;l++) //v'
				{
					Chiq[i][k*2*nv+l] += (Gk[j*2*nv + k] * Gk[k1*2*nv + l])/norm;
				}
			}
		}
	}
	
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

int Cq2Calc(double * const Cq2, dcomp ** const Chiq, const dcomp * const a, const int nk, const int nv)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	
	int i,k;
	dcomp blubb;
	
	for(i=0;i<ndistk;i++) //q
	{
		
		blubb = 0.;
		for(k=-nv;k<nv;k++) //v
		{
			blubb += a[k]*a[k]*Chiq[i][(2*nv + 1)*k] / (1. - a[k]*a[k]*Chiq[i][(2*nv + 1)*k]);
		}
		Cq2[i] = 1. / (1. + real(blubb) );
	}
	
	return(0);
}



int FqCalc(dcomp ** const Fup, dcomp ** const Fdown, const double * const Cq2, dcomp ** const Chiq, const dcomp * const a, const int nk, const int nv)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	
	int i,k,l;
	
	for(i=0;i<ndistk;i++) //q
	{
		for(k=-nv;k<nv;k++) //v
		{
			for(l=-nv;l<nv;l++) //v'
			{
				Fup[i][k*2*nv+l] = Cq2[i]* ( a[k]*a[k] / (1. - a[k]*a[k]*Chiq[i][k*2*nv+k]) ) * ( a[l]*a[l] / (1. - a[l]*a[l]*Chiq[i][l*2*nv+l]) ) ;
				Fdown[i][k*2*nv+l] = a[k]*a[l] / (1. - a[k]*a[l]*Chiq[i][k*2*nv+l]);
			}
		}
	}
	
	return(0);
}


int SigmaCalc(dcomp * const DualSig , dcomp ** const Fup, dcomp ** const Fdown, const double * const Cq2, dcomp ** const Chiq, const dcomp * const a, dcomp * const Gk, const int nk , const int nv, const int* const ktoq, const int* const ksym, const int* const qtok, const int* const * const kmap, const int* const * const ksum, const int* const * const kdif)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	const int nnk = nk*nk;
	
	const  double norm = nnk;
	int i,j,k,l;
	
	int k1,q1,k2;
	
	
	for(j=0;j<nnk;j++) //k
	{
		for(k=-nv;k<nv;k++) //v
		{
			DualSig[2*nv*j + k] = 0.;
		}
	}
	
// 	cout << "out" << endl;
	for(i=0;i<ndistk;i++) //k
	{
		k1 = qtok[i];
		
		for(j=0;j<nnk;j++) //q
		{
			q1 = ktoq[j];
			k2 = ksum[k1][j];
			
//   cout << "  " << q1 << "/" << ndistk << "  " << j << "/" << nnk << "  " << k1 << "/" << nnk << "    " << ksum[k1][j] << "/" << nnk << endl;
			for(k=-nv;k<nv;k++) //v
			{
				for(l=-nv;l<nv;l++) //v'
				{
					DualSig[2*nv*k1 + k] -= ( Fup[q1][2*nv*k + l] + a[k]*a[k]*a[l]*a[l]*Chiq[q1][2*nv*k + l]/2. ) * Gk[k2*2*nv + l]/norm;
				}
				DualSig[2*nv*k1 + k] += ( Fdown[q1][2*nv*k + k] + a[k]*a[k]*a[k]*a[k]*Chiq[q1][2*nv*k + k]/2. ) * Gk[k2*2*nv + k]/norm;
			}
			
		}
	}
	for(i=0;i<ndistk;i++) //k
	{
		k1 = qtok[i];
		
		for(j=1;j<8;j++) //symmetries
		{
			k2 = kmap[j][k1];
			
			for(k=-nv;k<nv;k++) //v
			{
				DualSig[2*nv*k2 + k] = DualSig[2*nv*k1 + k];
			}
			
		}
	}
	
	return(0);
}



class DFphParams
{
	public:
	int nk, nkin, nv, nvin;
	int nnk;
	int ndistk;
	
	//one-particle quantities
	dcomp* Sigmadual;
	dcomp* Sigmadualold;
	dcomp* Sigmacor;
	dcomp* Sigmadualhf;
	dcomp* G0dual; 
	dcomp* Gdual; 
	dcomp* Gloc;
	dcomp* Gdualloc;
	
	//k-grid quantities
	int ** ksum;
	int ** kdif;
	int ** kmap;
	int * ksym;
	int * ktoq;
	int * qtok;
	
	//Vertex quantities
	dcomp* a;
	
	dcomp ** Fup;
	dcomp ** Fdown;
	
	dcomp ** Chiq;
	double * Cq2;
	
	DFphParams( const int innk , const int innv , const int innkin , const int innvin )
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
		
		Sigmadualold = new dcomp [2*nv*nnk]; 
		Sigmadualold += nv;
		
		Sigmacor = new dcomp [2*nv*nnk]; 
		Sigmacor += nv;
		
		Sigmadualhf = new dcomp [2*nv];
		Sigmadualhf += nv;
		
		Gloc = new dcomp [2*nv]; 
		Gloc += nv;
		
		Gdualloc = new dcomp [2*nv]; 
		Gdualloc += nv;
		
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
		
		Sigmadualold -= nv;
		delete Sigmadualold;
		
		Sigmacor -= nv;
		delete Sigmacor;
		
		Sigmadualhf -= nv;
		delete Sigmadualhf;
		
		Gloc -= nv;
		delete Gloc;
		
		Gdualloc -= nv;
		delete Gdualloc;
		
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
// 		{
// 			int i,j;
// 			for (i = 0; i< nnk; i++)
// 			{
//               for (j = 0; j< 1; j++)
//               {
//                   cout << i << "/" << nnk << "    " << j << "/" << nnk << "    " << ksum[i][j] << "/" << nnk << endl;
//               }
// 			}
// 		}
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
		int i;
		
		a = new dcomp [2*nv];
		a += nv;
        
		Fup = new dcomp *[ndistk];
		Fdown = new dcomp *[ndistk];
		
		Chiq = new dcomp *[ndistk];
		Cq2 = new double [ndistk];
        
        for (i = 0; i< ndistk; i++)
		{
			*(Fup + i) = new dcomp [4*nv*nv];
			*(Fdown + i) = new dcomp [4*nv*nv];
			*(Chiq + i) = new dcomp [4*nv*nv];
			
			*(Fup + i) += (2*nv + 1)*nv;
			*(Fdown + i) += (2*nv + 1)*nv;
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
			Fdown[i] -= (2*nv + 1)*nv;
			Chiq[i] -= (2*nv + 1)*nv;
			
			delete Fup[i];
			delete Fdown[i];
			delete Chiq[i];
		}
		
		delete Fup;
		delete Fdown;
		
		delete Chiq;
		delete Cq2;
		
		return(0);
	};
	
	int ReadDMFT()
	{
		int i, j;
		
		dcomp* ain = new dcomp [2*nvin];
		
		readbin ("a" , ain , 2*nvin );
		
		ain += nvin;
		
		for(i=-nv; i<nv;  i++)
		{
			a[i] = ain[i];
		}
		
		ain -= nvin;
		
		delete ain;
		
		dcomp* Gin = new dcomp [2*nvin];
		
		readbin ("G1" , Gin , 2*nvin );
		
		Gin += nvin;
		for(i=-nv; i<nv;  i++)
		{
			Gloc[i] = Gin[i];
		}
		Gin -= nvin;
		
		delete Gin;
		
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
	
	int LadderDualSigma()
	{
		int i,j;
		for(i=0; i< nnk; i++)
		{
			for(j=-nv; j< nv; j++)
			{
				Sigmadual[i*2*nv + j] = Sigmadual[i*2*nv + j] / (1. -  Sigmadual[i*2*nv + j] * Gloc[j]);
			}
		}
		
		return(0);
	}
	
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
	
	int ReadSigCors()
	{
		readbin ("SigmaCork" , Sigmacor-nv , 2*nv*nk*nk );
		
		return(0);
	}
	
	void WriteSigCors()
	{
		writebin ("SigmaCork" , Sigmacor-nv , 2*nv*nnk );
	}
	
	int ReadDualSig()
	{
		readbin("DualSig", Sigmadual-nv, 2*nv*nnk);
		
		return(0);
	}
	
	void WriteDualSig()
	{
		writebin ("DualSig" , Sigmadual-nv , 2*nv*nnk );
	}
	
	void WriteGdual()
	{
		writebin ("Gdual" , Gdual-nv , 2*nv*nnk );
	}
	
	
	int LadderCalc()
	{
		ChiCalc(Chiq, Gdual, nk, nv, qtok, kmap, ksum, kdif);
		
		Cq2Calc(Cq2, Chiq, a, nk, nv);
		
		FqCalc(Fup, Fdown, Cq2, Chiq, a, nk, nv);
		
		return(0);
	}
	
	int SigCalc()
	{
//         cout << "   " << endl;
		SigmaCalc(Sigmadual , Fup, Fdown, Cq2, Chiq, a, Gdual, nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);		
		
		return(0);
	}
	
	//for ladder calculations in channel cha, 0 is ph, 1 is pp
	int LadderAndSigCalc()
	{
		//call Parquetiter with both ph and pp contributions to get right Sigma
		LadderCalc();
		SigCalc();
		
		//call Parquetiterph/pp to calculate actual vertex
		return(0);	
	}
};
