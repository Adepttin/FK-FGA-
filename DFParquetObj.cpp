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

//class for self contained DF parquet calculation based on DMFT results

//------------------------
//Bethe-Salpeter equations
//------------------------

//Bethe-Salpeter equation in the ph channel for fixed k and k'
int kBSph(dcomp * const Pphup, dcomp * const Pphdown, dcomp ** const Gphup, dcomp ** const Gphdown, dcomp ** const Fup, dcomp ** const Fdown , const int nk , const int nv)
{	
	const int nnk = nk*nk;
	const int nnv = 2*nv + 1;
	int i,j,k,l,m;
	int c1;
	int cv;
	
	//initialize Phis to zero
	for(i=-nv; i<nv;  i++)
	{
		for(j=-nv; j<nv;  j++)
		{
			Pphup[i*2*nv + j] = 0.;
			Pphdown[i*2*nv + j] = 0.;
		}
	}
	
	for(i=0; i<nk; i++)
	{
		for(j=0; j<nk; j++)
		{
			c1 = i*nk + j;
			
			for(k=-nv; k<nv; k++)
			{
				for(l=-nv; l<nv; l++)
				{					
					cv = k*2*nv + l;
					
					//calculate Phi down
					Pphdown[cv] -= (Fdown[c1][cv] * Gphdown[c1*nnk][cv]);
					
					//calculate Phi up
					Pphup[cv] += (Fup[c1][cv] * Gphdown[c1*nnk][l*nnv] + Fdown[c1][k*nnv] * Gphup[c1*nnk][cv]);
					
					for(m=-nv; m<nv; m++)
					{
						
						Pphup[cv] -= (Fup[c1][k*2*nv+m] * Gphup[c1*nnk][m*2*nv+l]);
						
					}
				}
			}
		}
	}
	
	//normalization 1/(nk)^2
	for(i=-nv; i<nv;  i++)
	{
		for(j=-nv; j<nv;  j++)
		{
			Pphup[i*2*nv + j] /= nnk;
			Pphdown[i*2*nv + j] /= nnk;
		}
	}
	
	return (0);
}

//Bethe-Salpeter equation in the ph channel with vertex functions at fixed q
int qBSph(dcomp ** const Pphup, dcomp ** const Pphdown, dcomp ** const Gphup, dcomp ** const Gphdown, dcomp ** const Fup, dcomp ** const Fdown , const int nk , const int nv)
{
	const int nnk = nk*nk;
	
	int i,j,k,l;
	int c1,c2;
	
	for(i=0; i<nk;  i++) //k
	{
		for(j=0; j<nk;  j++)
		{
			//k
			c1 = i*nk + j;

			for(k=0; k<nk;  k++) //k'
			{
				for(l=0; l<nk;  l++)
				{
					//k'
					c2 = k*nk + l;
					
					//call Bethe-Salpeter equation for fixed k and k'
					kBSph( Pphup[c1*nnk+c2], Pphdown[c1*nnk+c2], (Gphup+c2), (Gphdown+c2), (Fup+c1*nnk), (Fdown+c1*nnk), nk , nv);
				}
			}
		}
	}
	
	return (0);
}

//Bethe-Salpeter equation in the ph channel
int BSph(dcomp *** const Pphup, dcomp *** const Pphdown, dcomp *** const Gphup, dcomp *** const Gphdown, dcomp *** const Fup, dcomp *** const Fdown , const int nk , const int nv)
{	
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	int i;
	
	for(i=0; i<ndistk;  i++) //q
	{
		//call qBSph for fixed q
		qBSph( Pphup[i], Pphdown[i], Gphup[i], Gphdown[i], Fup[i], Fdown[i], nk , nv);
	}
	
	return (0);
}



//Bethe-Salpeter equation in the pp channel
int BSpp(dcomp *** const Pppup, dcomp *** const Pppdown, dcomp *** const Gppup, dcomp *** const Gppdown, dcomp *** const Fup, dcomp *** const Fdown , const int nk , const int nv, const int* const ktoq, const int* const ksym, const int* const qtok, const int* const * const kmap, const int* const * const ksum, const int* const * const kdif)
{	
	const int ndistk = (nk/2 + 1)*(nk/2+2)/2;
	const int nnk = nk*nk;
	
	int i,j,k,l,m,n;
	int kq,qkp,qkkp,q1,q2,khelp;
	int s1,s2;
	int dumk1, dumk2;
	
	for(i = 0; i < ndistk; i++) //q
	{	
		//q as k variable (on full BZ)
		kq = qtok[i];
		
		for(j = 0; j < nnk; j++)  //k
		{
			for(k = 0; k < nnk; k++)  //k'
			{
				//k + q
				qkp = ksum[kq][k];
				//k + k' + q
				qkkp = ksum[j][qkp];
				
				//initialize Phis to zero
				for(m = -nv; m < nv; m++)  //v
				{
					for(n = -nv; n < nv; n++)  //vp
					{
						Pppup[i][j*nnk+k][m*2*nv+n] = 0.;
						Pppdown[i][j*nnk+k][m*2*nv+n] = 0.;
					}
				}
				
				for(l = 0; l < nnk; l++)  //k1
				{
					//k + q - k1
					dumk1 = kdif[qkp][l];
					//k + q - k1 as q index (on irreducible BZ)
					q1 = ktoq[dumk1];
					//symmetry index for k + q - k1
					s1 = ksym[dumk1];
					//k1 - k'
					dumk2 = kdif[l][k];
					//k1 - k' as q index (on irreducible BZ)
					q2 = ktoq[dumk2];
					//symmetry index for k1 - k'
					s2 = ksym[dumk2];
					//k + k' + q - k1
					khelp = kdif[qkkp][l];
					
					//index for (k,k1)
					dumk1 = (kmap[s1][j] * nnk + kmap[s1][l]);
					//index for (k + k' + q -k1, k')
					dumk2 = (kmap[s2][khelp] * nnk + kmap[s2][k]);
					
					for(m = -nv; m < nv; m++)  //v
					{
						for(n = -nv; n < nv; n++)  //vp
						{
							//calculate Phi up
							Pppup[i][j*nnk+k][m*2*nv+n] += 0.5*( Gppup[q1][dumk1][m*2*nv+n] * Fup[q2][dumk2][m*2*nv+n] + Gppdown[q1][dumk1][m*2*nv+n] * Fdown[q2][dumk2][n*2*nv+m] );
							
							//calculate Phi down
							Pppdown[i][j*nnk+k][m*2*nv+n] += 0.5*( Gppup[q1][dumk1][m*2*nv+n] * Fdown[q2][dumk2][m*2*nv+n] + Gppdown[q1][dumk1][m*2*nv+n] * Fup[q2][dumk2][n*2*nv+m] );
						}
					}
					
				}
				
				//normalization 1/(nk)^2
				for(m = -nv; m < nv; m++)  //v
				{
					for(n = -nv; n < nv; n++)  //vp
					{
						Pppup[i][j*nnk+k][m*2*nv+n] /= (1.*nnk);
						Pppdown[i][j*nnk+k][m*2*nv+n] /= (1.*nnk);
					}
				}
				
			}
		}
		
	}
	
	return (0);
}

//------------------------
//parquet equation
//------------------------

//parquet equation
int ParquetReassembly(dcomp *** const Gppup, dcomp *** const Gppdown, dcomp *** const Gphup, dcomp *** const Gphdown, dcomp *** const Fup, dcomp *** const Fdown, dcomp *** const Pppup, dcomp *** const Pppdown, dcomp *** const Pphup, dcomp *** const Pphdown, dcomp * const Lambdaup, dcomp * const Lambdadown, dcomp * const Gk, const int nk , const int nv, const int* const ktoq, const int* const ksym, const int* const qtok, const int* const * const kmap, const int* const * const ksum, const int* const * const kdif)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	const int nnk = nk*nk;
	
	int i,j,k,l,m;
	int q1,k1;
	int s1;
	int dumk1, dumk2;
	
	for(i=0;i<ndistk;i++) //q
	{
		for(j=0;j<nnk;j++) //k
		{
			//k + q
			k1 = ksum[qtok[i]][j];
			
			for(k=0;k<nnk;k++) //k'
			{
				//k' - k
				dumk1 = (kdif[k][j]);
				//k' - k as q index on irreducible BZ
				q1 = ktoq[dumk1];
				//symmetry index for k' - k
				s1 = ksym[dumk1];
				
				//index for (k,k')
				dumk1 = j*nnk+k;
				//index for (k,k+q)
				dumk2 = (kmap[s1][j])*nnk+kmap[s1][k1];
				
				for(l=-nv;l<nv;l++) //v
				{
					for(m=-nv;m<nv;m++) //v'
					{
						//calculate F up
						Fup[i][dumk1][l*2*nv+m] = Lambdaup[l*2*nv+m] + Pppup[i][dumk1][l*2*nv+m] + Pphup[i][dumk1][l*2*nv+m] + Pphdown[q1][dumk2][l*2*nv+m];
						
						//calculate F down
						Fdown[i][dumk1][l*2*nv+m] = Lambdadown[l*2*nv+m] + Pppdown[i][dumk1][l*2*nv+m] + Pphdown[i][dumk1][l*2*nv+m] + Pphup[q1][dumk2][l*2*nv+m];
						
						//calculate Gamma pp up
						Gppup[i][dumk1][l*2*nv+m] = (Fup[i][dumk1][l*2*nv+m] - Pppup[i][dumk1][l*2*nv+m])*(Gk[k*2*nv + m])*(Gk[k1*2*nv + l]);
						
						//calculate Gamma pp down
						Gppdown[i][dumk1][l*2*nv+m] = (Fdown[i][dumk1][l*2*nv+m] - Pppdown[i][dumk1][l*2*nv+m])*(Gk[k*2*nv + l])*(Gk[k1*2*nv + m]);
						
						//calculate Gamma ph up
						Gphup[i][dumk1][l*2*nv+m] = (Fup[i][dumk1][l*2*nv+m] - Pphup[i][dumk1][l*2*nv+m])*(Gk[j*2*nv + m])*(Gk[k1*2*nv + m]);
						
						//calculate Gamma ph down
						Gphdown[i][dumk1][l*2*nv+m] = (Fdown[i][dumk1][l*2*nv+m] - Pphdown[i][dumk1][l*2*nv+m])*(Gk[j*2*nv + l])*(Gk[k1*2*nv + m]);
						
					}
				}
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

//prepare calculation of dual self energy
int InitSigSummand(dcomp*** const SigSummand, dcomp *** const Fup, dcomp * const Lambdaup, dcomp * const Gk, const int nk , const int nv, const int* const qtok, const int* const * const ksum)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	const int nnk = nk*nk;
	int i,j,k,l,m;
	int qeff;
	int kpq;
	
	dcomp* Leff = new dcomp [4*nv*nv];
	
	Leff += nv*(2*nv+1);
	
	//prepare reordered and dediagonalized Lambda (local F) including k-normalization here
	{
		
		for(i=-nv; i<nv;  i++) //v
		{
			for(j=-nv; j<nv; j++) //v'
			{
				Leff[i*2*nv + j] = Lambdaup[j*2*nv + i]/double(nnk);
			}
			
			Leff[i*2*nv + i] = 0.;
		}
	}
	
	//initialize SigSummand to zero
	for(i=0; i<ndistk; i++) //q
	{
		for(j=0; j<nnk; j++) //k
		{
			for(l=0; l<nv; l++) //v
			{
				SigSummand[i][j][l] = 0.;
			}
		}
	}
	
	for(i=0; i<ndistk; i++) //q
	{
		//q as k index on full BZ
		qeff = qtok[i];
		
		for(j=0; j<nnk; j++) //k
		{
			for(k=0; k<nnk; k++) //k'
			{
				//k' + q
				kpq = ksum[qeff][k];
				
				for(l=0; l<nv; l++) //v
				{
					for(m=-nv; m<nv; m++) //v'
					{
						//calculate SigSummand
						SigSummand[i][j][l] -= ( Fup[i][nnk*j + k][2*l*nv + m] )*Gk[2*nv*k + m]*Gk[2*nv*kpq + m]*( Leff[2*l*nv + m] );
					}
					
				}
				
			}
		}
	}
	
	//release memory
	Leff -= nv*(2*nv+1);
	delete Leff;
	
	return(0);
}

//calculate dual self energy via equation of motion
//without Hartree term
int CalcSigmaDual(dcomp* const Sigmadual, dcomp*** const SigSummand, dcomp * const Gk, const int nk , const int nv, const int* const ktoq, const int* const ksym, const int* const * const kmap, const int* const * const ksum)
{
	const int nnk = nk*nk;
	int i,j,l;
	int qeff,keff;
	int s;
	int kpq;
	
	//initialise dual self energy to zero
	for(j=0; j<nnk; j++) //k
	{
		for(l=0; l<nv; l++) //v
		{
			Sigmadual[2*nv*j + l] = 0.;
		}
	}
	
	for(i=0; i<nnk; i++) //k
	{
		//symmetry index of k
		s = ksym[i];
		//k as q index on irreducible BZ
		qeff = ktoq[i];
		
		for(j=0; j<nnk; j++) //k'
		{
			//k' as index on irreducible BZ
			keff = kmap[s][j];
			//k + k'
			kpq = ksum[i][j];
			
			for(l=0; l<nv; l++) //v
			{
				//calculation of dual self energy
				Sigmadual[2*nv*j + l] += SigSummand[qeff][keff][l] * Gk[2*nv*kpq + l];
			}
		}
	}
	
	
	for(j=0; j<nnk; j++)
	{
		for(l=0; l<nv; l++)
		{
			Sigmadual[2*nv*j + l] /= nnk;
			Sigmadual[2*nv*j - (l+1)] = conj(Sigmadual[2*nv*j + l]);
		}
	}

	return(0);
}

//calculate the Hartree Fock term (is purely local) for the self energy
int CalcSigmaDualHF(dcomp * const Gdual, dcomp * const Gdualloc, dcomp * const Flocup, dcomp * const Flocdown, dcomp * const Sigmadualhf, const int nk, const int nv)
{
	const int nnk = nk*nk;
	const double norm = nnk;    
	int i, j, k;
	
	//initialise local dual propagator and Hartree term to zero
	for(i=-nv; i<nv; i++) //v
	{
		Gdualloc[i] = 0.;
		Sigmadualhf[i] = 0.;
	}
	
	//calculate local dual propagator as sum over BZ
	for(k=0; k<nnk; k++) //k
	{
		for(i=-nv; i<nv; i++) //v
		{
			Gdualloc[i] += Gdual[2*nv*k + i]/norm;
		}
	}
	
	//calculate Hartree self energy
	for(i=-nv; i<nv; i++) //v
	{
		for(j=-nv; j<nv; j++) //v'
		{
			Sigmadualhf[i] -= Flocup[2 * i*nv + j] * Gdualloc[j];
		}
		
		Sigmadualhf[i] += Flocup[2 * i*nv + i] * Gdualloc[i];
	}
		
	return(0);	
}

class DFParquetParams
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
	//auxiliary variable for dual self energy
	dcomp* Sigmadualold;
	//resulting DF corrections Sigma_corr(k,v)
	dcomp* Sigmacor;
	//Hartree term of dual self energy Sigma_HF(v)
	dcomp* Sigmadualhf;
	//bare dual propagator G_0(k,v)
	dcomp* G0dual; 
	//dual propagator G(k,v)
	dcomp* Gdual; 
	//local propagator of real fermions G_loc(v)
	dcomp* Gloc;
	//local propagator of dual fermions G_loc(v)
	dcomp* Gdualloc;
	
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
	
	//vertex quantities
	
	//local vertex function F_loc(v,v')
	dcomp* Flocup;
	dcomp* Flocdown;
	
	//full vertex function F(k,k',q)(v,v')
	dcomp *** Fup;
	dcomp *** Fdown;
	
	//irreducible vertex functions Gamma_ph(k,k',q)(v,v')
	dcomp *** Gphup;
	dcomp *** Gphdown;
	
	//reducible vertex functions Phi_ph(k,k',q)(v,v')
	dcomp *** Pphup;
	dcomp *** Pphdown;
	
	//irreducible vertex functions Gamma_pp(k,k',q)(v,v')
	dcomp *** Gppup;
	dcomp *** Gppdown;
	
	//reducible vertex functions Phi_pp(k,k',q)(v,v')
	dcomp *** Pppup;
	dcomp *** Pppdown;
	
	//little helper for calculation of dual self energy
	dcomp*** SigSummand;
	
	//constructor
	DFParquetParams( const int innk , const int innv , const int innkin , const int innvin )
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
		
		Flocup = new dcomp [4*nv*nv];
		Flocdown = new dcomp [4*nv*nv];
		
		Flocup += nv*(2*nv+1);
		Flocdown += nv*(2*nv+1);
		
		Fup = new dcomp **[ndistk];
		Fdown = new dcomp **[ndistk];
		
		Gphup = new dcomp **[ndistk];
		Gphdown = new dcomp **[ndistk];
		
		Pphup = new dcomp **[ndistk];
		Pphdown = new dcomp **[ndistk];
		
		
		Gppup = new dcomp **[ndistk];
		Gppdown = new dcomp **[ndistk];
		
		Pppup = new dcomp **[ndistk];
		Pppdown = new dcomp **[ndistk];
		
		SigSummand = new dcomp **[ndistk];
		
		for (i = 0; i< ndistk; i++)
		{
			*(Fup + i) = new dcomp *[nnk*nnk];
			*(Fdown + i) = new dcomp *[nnk*nnk];
			
			
			*(Gphup + i) = new dcomp *[nnk*nnk];
			*(Gphdown + i) = new dcomp *[nnk*nnk];
			
			*(Pphup + i) = new dcomp *[nnk*nnk];
			*(Pphdown + i) = new dcomp *[nnk*nnk];
			
			
			*(Gppup + i) = new dcomp *[nnk*nnk];
			*(Gppdown + i) = new dcomp *[nnk*nnk];
			
			*(Pppup + i) = new dcomp *[nnk*nnk];
			*(Pppdown + i) = new dcomp *[nnk*nnk];
			
			*(SigSummand + i) = new dcomp *[nnk];
			
			for (j = 0; j < nnk*nnk; j++)
			{
				Fup[i][j] = new dcomp [4*nv*nv];
				Fdown[i][j] = new dcomp [4*nv*nv];
				
				Fup[i][j] += (2*nv + 1)*nv;
				Fdown[i][j] += (2*nv + 1)*nv;
				
				Gphup[i][j] = new dcomp [4*nv*nv];
				Gphdown[i][j] = new dcomp [4*nv*nv];
				
				Gphup[i][j] += (2*nv + 1)*nv;
				Gphdown[i][j] += (2*nv + 1)*nv;
				
				Pphup[i][j] = new dcomp [4*nv*nv];
				Pphdown[i][j] = new dcomp [4*nv*nv];
				
				Pphup[i][j] += (2*nv + 1)*nv;
				Pphdown[i][j] += (2*nv + 1)*nv;
				
				
				Gppup[i][j] = new dcomp [4*nv*nv];
				Gppdown[i][j] = new dcomp [4*nv*nv];
				
				Gppup[i][j] += (2*nv + 1)*nv;
				Gppdown[i][j] += (2*nv + 1)*nv;
				
				Pppup[i][j] = new dcomp [4*nv*nv];
				Pppdown[i][j] = new dcomp [4*nv*nv];
				
				Pppup[i][j] += (2*nv + 1)*nv;
				Pppdown[i][j] += (2*nv + 1)*nv;
			}
			
			for (j = 0; j < nnk; j++)
			{
				SigSummand[i][j] = new dcomp [nv];
			}
		}
		
		return(0);
	};
	
	int DeleteVertexStorage()
	{
		int i, j;
		
		for (i = 0; i< ndistk; i++)
		{
			
			for (j = 0; j < nnk*nnk; j++)
			{
				Fup[i][j] -= (2*nv + 1)*nv;
				Fdown[i][j] -= (2*nv + 1)*nv;
				delete Fup[i][j];
				delete Fdown[i][j];
				
				
				Gphup[i][j] -= (2*nv + 1)*nv;
				Gphdown[i][j] -= (2*nv + 1)*nv;
				delete Gphup[i][j];
				delete Gphdown[i][j];
				
				Pphup[i][j] -= (2*nv + 1)*nv;
				Pphdown[i][j] -= (2*nv + 1)*nv;
				delete Pphup[i][j];
				delete Pphdown[i][j];
				
				
				Gppup[i][j] -= (2*nv + 1)*nv;
				Gppdown[i][j] -= (2*nv + 1)*nv;
				delete Gppup[i][j];
				delete Gppdown[i][j];
				
				Pppup[i][j] -= (2*nv + 1)*nv;
				Pppdown[i][j] -= (2*nv + 1)*nv;
				delete Pppup[i][j];
				delete Pppdown[i][j];
			}
			for (j = 0; j < nnk; j++)
			{
				delete SigSummand[i][j];
			}
			
			delete Fup[i];
			delete Fdown[i];
			
			
			delete Gphup[i];
			delete Gphdown[i];
			
			delete Pphup[i];
			delete Pphdown[i];
			
			
			delete Gppup[i];
			delete Gppdown[i];
			
			delete Pppup[i];
			delete Pppdown[i];
			
			delete SigSummand[i];
		}
		
		delete Fup;
		delete Fdown;
		
		
		delete Gphup;
		delete Gphdown;
		
		delete Pphup;
		delete Pphdown;
		
		
		delete Gppup;
		delete Gppdown;
		
		delete Pppup;
		delete Pppdown;
		
		delete SigSummand;
		
		return(0);
	};
	
	//read in DMFT results
	int ReadDMFT()
	{
		//read in local vertex function F
		
		int i, j;
		
		dcomp* Fupin = new dcomp [4*nvin*nvin];
		dcomp* Fdownin = new dcomp [4*nvin*nvin];
		
		readbin ("LambdaUp" , Fupin , 4*nvin*nvin );
		readbin ("LambdaDown" , Fdownin , 4*nvin*nvin );
		
		Fupin += nvin*(2*nvin+1);
		Fdownin += nvin*(2*nvin+1);
		
		for(i=-nv; i<nv;  i++)
		{
			for(j=-nv; j<nv; j++)
			{
				Flocup[i*2*nv + j] = Fupin[i*2*nvin + j];
				Flocdown[i*2*nv + j] = Fdownin[i*2*nvin + j];
			}
		}
		
		Fupin -= nvin*(2*nvin+1);
		Fdownin -= nvin*(2*nvin+1);
		
		delete Fupin;
		delete Fdownin;
		
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
	
	//reset vertex functions
	int ResetVertex()
	{
		int i,j,k;
		
		//initialise reducible Phis to zero
		for(i=0; i<ndistk; i++)
		{
			for(j=0; j<nnk*nnk; j++)
			{
				for(k=-nv*(2*nv + 1); k<(nv-1)*(2*nv + 1)+1; k++)
				{
					
					Pphup[i][j][k] = 0.;
					Pphdown[i][j][k] = 0.;
					
					Pppup[i][j][k] = 0.;
					Pppdown[i][j][k] = 0.;
					
				}
			}
		}
		
		//call the parquet equation to initialise F and Gammas to local F
		this->Parquetiter();
		
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
				//Dyson equation
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
	
	//call Bethe-Salpeter equations
	int BSiter()
	{		
		int retval = 0;
		
		//for ph channel
		retval += BSph( Pphup, Pphdown, Gphup, Gphdown, Fup, Fdown , nk , nv);
		//for pp channel
		retval += BSpp( Pppup, Pppdown, Gppup, Gppdown, Fup, Fdown , nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);
		
		return (retval);
	}
	
	//call Bethe-Salepter equation only for ph channel
	int BSiterph()
	{
		int retval = 0;
		
		retval += BSph( Pphup, Pphdown, Gphup, Gphdown, Fup, Fdown , nk , nv);
		
		return (retval);
	}
	
	//call Bethe-Salepter equation only for pp channel
	int BSiterpp()
	{
		int retval = 0;
		
		retval += BSpp( Pppup, Pppdown, Gppup, Gppdown, Fup, Fdown , nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);
		
		return (retval);
	}
	
	//call parquet equation
	int Parquetiter()
	{
		ParquetReassembly(Gppup, Gppdown, Gphup, Gphdown, Fup, Fdown, Pppup, Pppdown, Pphup, Pphdown, Flocup, Flocdown, Gdual, nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);
		return(0);
	}
	
	//calculate dual self energy via equation of motion
	int SigCalc()
	{
		
		int i,k;
		
		//initialise the auxiliary quantity SigSummand
		InitSigSummand(SigSummand, Fup, Flocup, Gdual, nk , nv, qtok, ksum);	
		dcomp* dummy;
		
		//store old dual self energy
		dummy = Sigmadualold;
		Sigmadualold = Sigmadual;
		Sigmadual = dummy;
		
		//calculate dual self energy without Hartree term
		CalcSigmaDual(Sigmadual, SigSummand, Gdual, nk , nv, ktoq, ksym, kmap, ksum);
		
		//calculate Hartree term
		CalcSigmaDualHF(Gdual, Gdualloc, Flocup, Flocdown, Sigmadualhf, nk, nv);
		
		//write Hartree term
		//writebin("Sigmadualhf", Sigmadualhf-nv, 2*nv);
		
		//write local dual propagator
		//writebin("Gdualloc", Gdualloc-nv, 2*nv);
		
		//add Hartree term to dual self energy
		for(k=0; k<nnk; k++) //k
		{
			for(i=-nv; i<nv; i++) //v
			{
				Sigmadual[2*nv*k+i]=Sigmadualhf[i]+Sigmadual[2*nv*k+i];
			}
		}
		
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
