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

// Self contained DF based on DMFT results


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

int kBSph(dcomp * const Pphup, dcomp * const Pphdown, dcomp ** const Gphup, dcomp ** const Gphdown, dcomp ** const Fup, dcomp ** const Fdown , const int nk , const int nv)
{	
	const int nnk = nk*nk;
	const int nnv = 2*nv + 1;
	int i,j,k,l,m;
	int c1;
	int cv;
	
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
					Pphdown[cv] -= (Fdown[c1][cv] * Gphdown[c1*nnk][cv]);
					
					Pphup[cv] += (Fup[c1][cv] * Gphdown[c1*nnk][l*nnv] + Fdown[c1][k*nnv] * Gphup[c1*nnk][cv]);
					
					for(m=-nv; m<nv; m++)
					{
						
						Pphup[cv] -= (Fup[c1][k*2*nv+m] * Gphup[c1*nnk][m*2*nv+l]);
						
					}
				}
			}
		}
	}
	
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

int qBSph(dcomp ** const Pphup, dcomp ** const Pphdown, dcomp ** const Gphup, dcomp ** const Gphdown, dcomp ** const Fup, dcomp ** const Fdown , const int nk , const int nv)
{
	const int nnk = nk*nk;
	
	int i,j,k,l;
	int c1,c2;
	
	for(i=0; i<nk;  i++)
	{
		for(j=0; j<nk;  j++)
		{
			c1 = i*nk + j;
			for(k=0; k<nk;  k++)
			{
				for(l=0; l<nk;  l++)
				{
					
					c2 = k*nk + l;
					kBSph( Pphup[c1*nnk+c2], Pphdown[c1*nnk+c2], (Gphup+c2), (Gphdown+c2), (Fup+c1*nnk), (Fdown+c1*nnk), nk , nv);
				}
			}
		}
	}
	
	return (0);
}

int BSph(dcomp *** const Pphup, dcomp *** const Pphdown, dcomp *** const Gphup, dcomp *** const Gphdown, dcomp *** const Fup, dcomp *** const Fdown , const int nk , const int nv)
{	
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	int i;
	
	for(i=0; i<ndistk;  i++)
	{
		qBSph( Pphup[i], Pphdown[i], Gphup[i], Gphdown[i], Fup[i], Fdown[i], nk , nv);
	}
	
	return (0);
}


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
		kq = qtok[i];
		for(j = 0; j < nnk; j++)  //k
		{
			for(k = 0; k < nnk; k++)  //k'
			{
				qkp = ksum[kq][k];
				qkkp = ksum[j][qkp];
				
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
					dumk1 = kdif[qkp][l];
					q1 = ktoq[dumk1];
					s1 = ksym[dumk1];
					dumk2 = kdif[l][k];
					q2 = ktoq[dumk2];
					s2 = ksym[dumk2];
					khelp = kdif[qkkp][l];
					
					dumk1 = (kmap[s1][j] * nnk + kmap[s1][l]);
					dumk2 = (kmap[s2][khelp] * nnk + kmap[s2][k]);
					
					for(m = -nv; m < nv; m++)  //v
					{
						for(n = -nv; n < nv; n++)  //vp
						{
							Pppup[i][j*nnk+k][m*2*nv+n] += 0.5*( Gppup[q1][dumk1][m*2*nv+n] * Fup[q2][dumk2][m*2*nv+n] + Gppdown[q1][dumk1][m*2*nv+n] * Fdown[q2][dumk2][n*2*nv+m] );
							
							Pppdown[i][j*nnk+k][m*2*nv+n] += 0.5*( Gppup[q1][dumk1][m*2*nv+n] * Fdown[q2][dumk2][m*2*nv+n] + Gppdown[q1][dumk1][m*2*nv+n] * Fup[q2][dumk2][n*2*nv+m] );
						}
					}
					
				}
				
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
			k1 = ksum[qtok[i]][j];
			
			for(k=0;k<nnk;k++) //k'
			{
				dumk1 = (kdif[k][j]);
				q1 = ktoq[dumk1];
				s1 = ksym[dumk1];
				
				dumk1 = j*nnk+k;
				dumk2 = (kmap[s1][j])*nnk+kmap[s1][k1];
				
				for(l=-nv;l<nv;l++) //v
				{
					for(m=-nv;m<nv;m++) //v'
					{
						Fup[i][dumk1][l*2*nv+m] = Lambdaup[l*2*nv+m] + Pppup[i][dumk1][l*2*nv+m] + Pphup[i][dumk1][l*2*nv+m] + Pphdown[q1][dumk2][l*2*nv+m];
						
						Fdown[i][dumk1][l*2*nv+m] = Lambdadown[l*2*nv+m] + Pppdown[i][dumk1][l*2*nv+m] + Pphdown[i][dumk1][l*2*nv+m] + Pphup[q1][dumk2][l*2*nv+m];
						
						
						Gppup[i][dumk1][l*2*nv+m] = (Fup[i][dumk1][l*2*nv+m] - Pppup[i][dumk1][l*2*nv+m])*(Gk[k*2*nv + m])*(Gk[k1*2*nv + l]);
						
						Gppdown[i][dumk1][l*2*nv+m] = (Fdown[i][dumk1][l*2*nv+m] - Pppdown[i][dumk1][l*2*nv+m])*(Gk[k*2*nv + l])*(Gk[k1*2*nv + m]);
						
						
						Gphup[i][dumk1][l*2*nv+m] = (Fup[i][dumk1][l*2*nv+m] - Pphup[i][dumk1][l*2*nv+m])*(Gk[j*2*nv + m])*(Gk[k1*2*nv + m]);
						
						Gphdown[i][dumk1][l*2*nv+m] = (Fdown[i][dumk1][l*2*nv+m] - Pphdown[i][dumk1][l*2*nv+m])*(Gk[j*2*nv + l])*(Gk[k1*2*nv + m]);
						
					}
				}
			}
		}
	}
	
	return(0);
}

//calculates only contributions from the ph channel
int ParquetReassemblyph(dcomp *** const Gppup, dcomp *** const Gppdown, dcomp *** const Gphup, dcomp *** const Gphdown, dcomp *** const Fup, dcomp *** const Fdown, dcomp *** const Pppup, dcomp *** const Pppdown, dcomp *** const Pphup, dcomp *** const Pphdown, dcomp * const Lambdaup, dcomp * const Lambdadown, dcomp * const Gk, const int nk , const int nv, const int* const ktoq, const int* const ksym, const int* const qtok, const int* const * const kmap, const int* const * const ksum, const int* const * const kdif)
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
			k1 = ksum[qtok[i]][j];
			
			for(k=0;k<nnk;k++) //k'
			{
				dumk1 = (kdif[k][j]);
				q1 = ktoq[dumk1];
				s1 = ksym[dumk1];
				
				dumk1 = j*nnk+k;
				dumk2 = (kmap[s1][j])*nnk+kmap[s1][k1];
				
				for(l=-nv;l<nv;l++) //v
				{
					for(m=-nv;m<nv;m++) //v'
					{
						Fup[i][dumk1][l*2*nv+m] = Lambdaup[l*2*nv+m] + Pphup[i][dumk1][l*2*nv+m];
						Fdown[i][dumk1][l*2*nv+m] = Lambdadown[l*2*nv+m] + Pphdown[i][dumk1][l*2*nv+m];					
						
						Gphup[i][dumk1][l*2*nv+m] = (Fup[i][dumk1][l*2*nv+m] - Pphup[i][dumk1][l*2*nv+m])*(Gk[j*2*nv + m])*(Gk[k1*2*nv + m]);
						Gphdown[i][dumk1][l*2*nv+m] = (Fdown[i][dumk1][l*2*nv+m] - Pphdown[i][dumk1][l*2*nv+m])*(Gk[j*2*nv + l])*(Gk[k1*2*nv + m]);
						
					}
				}
			}
		}
	}
	
	return(0);
}

//calculates only contributions from the pp channel
int ParquetReassemblypp(dcomp *** const Gppup, dcomp *** const Gppdown, dcomp *** const Gphup, dcomp *** const Gphdown, dcomp *** const Fup, dcomp *** const Fdown, dcomp *** const Pppup, dcomp *** const Pppdown, dcomp *** const Pphup, dcomp *** const Pphdown, dcomp * const Lambdaup, dcomp * const Lambdadown, dcomp * const Gk, const int nk , const int nv, const int* const ktoq, const int* const ksym, const int* const qtok, const int* const * const kmap, const int* const * const ksum, const int* const * const kdif)
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
			k1 = ksum[qtok[i]][j];
			
			for(k=0;k<nnk;k++) //k'
			{
				dumk1 = (kdif[k][j]);
				q1 = ktoq[dumk1];
				s1 = ksym[dumk1];
				
				dumk1 = j*nnk+k;
				dumk2 = (kmap[s1][j])*nnk+kmap[s1][k1];
				
				for(l=-nv;l<nv;l++) //v
				{
					for(m=-nv;m<nv;m++) //v'
					{
						//Parquet equation only for channel pp
						Fup[i][dumk1][l*2*nv+m] = Lambdaup[l*2*nv+m] + Pppup[i][dumk1][l*2*nv+m];
						Fdown[i][dumk1][l*2*nv+m] = Lambdadown[l*2*nv+m] + Pppdown[i][dumk1][l*2*nv+m];
						
						//calculation only of Gamma_pp
						Gppup[i][dumk1][l*2*nv+m] = (Fup[i][dumk1][l*2*nv+m] - Pppup[i][dumk1][l*2*nv+m])*(Gk[k*2*nv + m])*(Gk[k1*2*nv + l]);					
						Gppdown[i][dumk1][l*2*nv+m] = (Fdown[i][dumk1][l*2*nv+m] - Pppdown[i][dumk1][l*2*nv+m])*(Gk[k*2*nv + l])*(Gk[k1*2*nv + m]);

						
					}
				}
			}
		}
	}
	
	return(0);
}


int InitSigSummand(dcomp*** const SigSummand, dcomp *** const Fup, dcomp * const Lambdaup, dcomp * const Gk, const int nk , const int nv, const int* const qtok, const int* const * const ksum)
{
	const int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	const int nnk = nk*nk;
	int i,j,k,l,m;
	int qeff;
	int kpq;
	
	dcomp* Leff = new dcomp [4*nv*nv];
	
	Leff += nv*(2*nv+1);
	
	//Prepare reordered and dediagonalized Lambda (k-normalization here)
	{
		
		for(i=-nv; i<nv;  i++)
		{
			for(j=-nv; j<nv; j++)
			{
				Leff[i*2*nv + j] = Lambdaup[j*2*nv + i]/double(nnk);
			}
			
			Leff[i*2*nv + i] = 0.;
		}
	}
	
	for(i=0; i<ndistk; i++)
	{
		for(j=0; j<nnk; j++)
		{
			for(l=0; l<nv; l++)
			{
				SigSummand[i][j][l] = 0.;
			}
		}
	}
	
	for(i=0; i<ndistk; i++)
	{
		qeff = qtok[i];
		
		for(j=0; j<nnk; j++)
		{
			for(k=0; k<nnk; k++)
			{
				kpq = ksum[qeff][k];
				
				for(l=0; l<nv; l++)
				{
					for(m=-nv; m<nv; m++)
					{
						SigSummand[i][j][l] -= ( Fup[i][nnk*j + k][2*l*nv + m] )*Gk[2*nv*k + m]*Gk[2*nv*kpq + m]*( Leff[2*l*nv + m] );
					}
					
				}
				
			}
		}
	}
	
	Leff -= nv*(2*nv+1);
	
	delete Leff;
	
	return(0);
}


int CalcSigmaDual(dcomp* const Sigmadual, dcomp*** const SigSummand, dcomp * const Gk, const int nk , const int nv, const int* const ktoq, const int* const ksym, const int* const * const kmap, const int* const * const ksum)
{
	const int nnk = nk*nk;
	int i,j,l;
	int qeff,keff;
	int s;
	int kpq;
	
	for(j=0; j<nnk; j++)
	{
		for(l=0; l<nv; l++)
		{
			Sigmadual[2*nv*j + l] = 0.;
		}
	}
	
	for(i=0; i<nnk; i++)
	{
		s = ksym[i];
		qeff = ktoq[i];
		
		for(j=0; j<nnk; j++)
		{
			keff = kmap[s][j];
			kpq = ksum[i][j];
			
			for(l=0; l<nv; l++)
			{
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

//calculate the Hartree Fock term for the self energy
int CalcSigmaDualHF(dcomp * const Gdual, dcomp * const Gdualloc, dcomp * const Flocup, dcomp * const Flocdown, dcomp * const Sigmadualhf, const int nk, const int nv)
{
	const int nnk = nk*nk;
	const double norm = nnk;    
	int i, j, k;
	
	for(i=-nv; i<nv; i++)
	{
		Gdualloc[i] = 0.;
		Sigmadualhf[i] = 0.;
	}
	
	//calculate local Gdual
	for(k=0; k<nnk; k++)
	{
		for(i=-nv; i<nv; i++)
		{
			Gdualloc[i] += Gdual[2 * nv*k + i]/norm;
		}
	}
	
	//calculate Hartree Fock self energy
	for(i=-nv; i<nv; i++)
	{
		for(j=-nv; j<nv; j++)
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
	dcomp* Flocup;
	dcomp* Flocdown;
	
	dcomp *** Fup;
	dcomp *** Fdown;
	
	
	dcomp *** Gphup;
	dcomp *** Gphdown;
	
	dcomp *** Pphup;
	dcomp *** Pphdown;
	
	
	dcomp *** Gppup;
	dcomp *** Gppdown;
	
	dcomp *** Pppup;
	dcomp *** Pppdown;
	
	//little helper
	dcomp*** SigSummand;
	
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
	
	int ReadDMFT()
	{
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
	
	int ResetVertex()
	{
		int i,j,k;
		
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
		
		this->Parquetiter();
		
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
	
	int BSiter()
	{		
		int retval = 0;
		
		retval += BSph( Pphup, Pphdown, Gphup, Gphdown, Fup, Fdown , nk , nv);
		retval += BSpp( Pppup, Pppdown, Gppup, Gppdown, Fup, Fdown , nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);
		
		return (retval);
	}
	
	int BSiterph()
	{
		int retval = 0;
		
		retval += BSph( Pphup, Pphdown, Gphup, Gphdown, Fup, Fdown , nk , nv);
		
		return (retval);
	}
	
	int BSiterpp()
	{
		int retval = 0;
		
		retval += BSpp( Pppup, Pppdown, Gppup, Gppdown, Fup, Fdown , nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);
		
		return (retval);
	}
	
	
	int Parquetiter()
	{
		
		ParquetReassembly(Gppup, Gppdown, Gphup, Gphdown, Fup, Fdown, Pppup, Pppdown, Pphup, Pphdown, Flocup, Flocdown, Gdual, nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);
		return(0);
	}
	
	int Parquetiterph()
	{
		
		ParquetReassemblyph(Gppup, Gppdown, Gphup, Gphdown, Fup, Fdown, Pppup, Pppdown, Pphup, Pphdown, Flocup, Flocdown, Gdual, nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);
		return(0);
	}
	
	int Parquetiterpp()
	{
		
		ParquetReassemblypp(Gppup, Gppdown, Gphup, Gphdown, Fup, Fdown, Pppup, Pppdown, Pphup, Pphdown, Flocup, Flocdown, Gdual, nk , nv, ktoq, ksym, qtok, kmap, ksum, kdif);
		return(0);
	}
	
	int SigCalc()
	{
		
		int i,k;
		
		InitSigSummand(SigSummand, Fup, Flocup, Gdual, nk , nv, qtok, ksum);	
		dcomp* dummy;
		
		dummy = Sigmadualold;
		Sigmadualold = Sigmadual;
		Sigmadual = dummy;
		
		CalcSigmaDual(Sigmadual, SigSummand, Gdual, nk , nv, ktoq, ksym, kmap, ksum);
		
		//add Hartree-Fock term
		CalcSigmaDualHF(Gdual, Gdualloc, Flocup, Flocdown, Sigmadualhf, nk, nv);
		
		writebin("Sigmadualhf", Sigmadualhf-nv, 2*nv);
		
		writebin("Gdualloc", Gdualloc-nv, 2*nv);
			
		for(k=0; k<nnk; k++)
		{
			for(i=-nv; i<nv; i++)
			{
				Sigmadual[2*nv*k+i]=Sigmadualhf[i]+Sigmadual[2*nv*k+i];
			}
		}
		
		return(0);
	}
	
	//for ladder calculations in channel cha, 0 is ph, 1 is pp
	int ParquetAndSigCalc(int cha)
	{
		//call Parquetiter with both ph and pp contributions to get right Sigma
		Parquetiter();
		SigCalc();
		
		//call Parquetiterph/pp to calculate actual vertex
		if (cha == 0)
		{
			Parquetiterph();
		}
		else if (cha == 1)
		{
			Parquetiterpp();
		}
		
		return(0);	
	}
};
