//calculate bubble term for chi
int calcBubbleChi(dcomp* const Gk, const double beta, const int nv, const int nk, int* const qtok, int** const ksum, dcomp* chibubble)
{
	int indexA, indexB;
	dcomp Gsum;
	int nnk = nk*nk;
	int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	
	for (int q=0; q < ndistk; q++)
	{
		for (int w=-nv; w < nv+1 ; w++)
		{
			
			Gsum = 0;
			
			//sum over k and v
			for (int k=0; k < nnk; k++)
			{
					for (int v=-nv; v < nv; v++)
					{
						//index for k and v
						indexA = 2*nv*k+v;
						
						//index for k+q and v+w
						indexB = 2*nv*ksum[qtok[q]][k]+(v+w);
						
						//account for v+w exceeding frequency range
						if( ((v+w) < nv) && ((v+w) >= -nv))
						{
							Gsum += Gk[indexA]*Gk[indexB];
						}
					}
			}
			
			chibubble[(2*nv+1)*q+w]=-Gsum/(double(nnk*nnk)*beta);
		}	
	}

	
	return(0);
}

//calculates vertex contribution to magnetic susceptibilty
//if FlocSep is not equal zero, then the first order diagram in Floc will be calculated separately on a larger frequency grid (nvFloc)
int calcVertexChi(dcomp* const Gk, dcomp*** const Fup, dcomp*** const Fdown,  dcomp * const Flocup,  dcomp * const Flocdown, const int FlocSep, const  double beta, const int nv, const int nk, const int nvFloc, int* const qtok, int** const ksum, dcomp* chivertex)
{
	int indexA, indexB, indexC, indexD;
	dcomp Fdownsum, Fupsum;
	int nnk = nk*nk;
	int ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	
	//calculating the first diagram (only Floc) separately on a nvFloc Matsubara frequency grid
	if (FlocSep)
	{
			
	}
	
	for (int q=0; q < ndistk; q++)
	{
		//calculating second term with Fdown
		for (int w=-nv; w < nv +1 ; w++)
		{
			Fdownsum=0;
			
			//summing up over k,k' and v
			for (int k=0; k < nnk; k++)
			{
				for (int kk=0; kk < nnk; kk++)
				{
					for (int v=-nv; v<nv; v++)
					{						
						//index for k and v
						indexA = 2*nv*k+v;
						
						//index for k+q and v+w
						indexB = 2*nv*ksum[qtok[q]][k]+(v+w);
						
						//index for k' and v
						indexC = 2*nv*kk+v;
						
						//index for k'+q and v+w
						indexD = 2*nv*ksum[qtok[q]][kk]+(v+w);
						
						
						//account for v+w exceeding frequency range
						if( ((v+w) < nv) && ((v+w) >= -nv))
						{
							Fdownsum -= Gk[indexA]*Gk[indexB]*Fdown[q][nnk*k+kk][2*nv*v+(v+w)]*Gk[indexC]*Gk[indexD];
						}
					}						
				}
			}
			
			chivertex[(2*nv+1)*q+w] = Fdownsum/(double(nnk*nnk)*beta*beta);
			
		}	
		
		//calculating first term with Fup (only for w=0)
		
		Fupsum = 0;
		
		//summing up over k,k' and v,v'
		for (int k=0; k < nnk; k++)
		{
			for (int kk=0; kk < nnk; kk++)
			{
				for (int v=-nv; v<nv; v++)
				{
					for (int vv=-nv; vv<nv; vv++)
					{
						//index for k and v
						indexA = 2*nv*k+v;

						//index for k+q and v
						indexB = 2*nv*ksum[qtok[q]][k]+v;

						//index for k' and v'
						indexC = 2*nv*kk+vv;

						//index for k'+q and v'
						indexD = 2*nv*ksum[qtok[q]][kk]+vv;

						Fupsum += Gk[indexA]*Gk[indexB]*Fup[q][nnk*k+kk][2*nv*v+vv]*Gk[indexC]*Gk[indexD];
					}
				}						
			}
		}
		
		chivertex[(2*nv+1)*q] += Fupsum/(double(nnk*nnk)*beta*beta);
	}
	
	return(0);
}


class MagnetismObject
{
	public:
	//grid quantities
	int nk, nv, nvin;
	int ndistk;
	
	//parameters
	double beta, mu;
	
	//dispersion
	double* Ek;
	
	//new index for sum of two k-vectors
	int ** ksum;
	//new index for q (IBZ) to k (full BZ)
	int * qtok;
	
	//Sigma from DMFT
	dcomp* Sigmaloc;
	//Sigma corrections from DF
	dcomp* Sigmacor;
	//full Green's function
	dcomp* Gk;
	
	//Matsubara frequencies
	dcomp* Matsus;
	
	//vertexfunctions
	dcomp* Flocup;
	dcomp* Flocdown;
	dcomp *** Fup;
	dcomp *** Fdown;
	
	//susceptibilites
	dcomp* chibubble;
	dcomp* chivertex;
	
	template <class Parquet>
	MagnetismObject(Parquet MyParquet, double inbeta, double inmu)
	{
		nk = MyParquet.nk;
		nv = MyParquet.nv;
		nvin = MyParquet.nvin;
		
		ksum = MyParquet.ksum;
		qtok = MyParquet.qtok;
		
		Sigmacor = MyParquet.Sigmacor;
		
		Flocup = MyParquet.Flocup;
		Flocdown = MyParquet.Flocdown;
		
		Fup = MyParquet.Fup;
		Fdown = MyParquet.Fdown;
		
		beta = inbeta;
		mu = inmu;
		
		ndistk = (nk/2 + 1)*(nk/2 + 2)/2;
	}
		
	int InitialiseStorage()
	{
		Ek = new double[nk*nk];
				
		Sigmaloc = new dcomp[2*nvin];
		Sigmaloc = Sigmaloc + nvin;
		
		Gk = new dcomp[2*nv*nk*nk];
		Gk = Gk + nv;
		
		chibubble = new dcomp [(2*nv+1)*ndistk];
		chibubble = chibubble + nv;
		
		chivertex = new dcomp [(2*nv+1)*ndistk];
		chivertex = chivertex + nv;
		
		return(0);
	}
	
	int DeleteStorage()
	{
		delete Ek;
		
		Sigmaloc = Sigmaloc - nvin;
		delete Sigmaloc;
		
		Gk = Gk - nv;
		delete Gk;
		
		chibubble = chibubble - nv;
		delete chibubble;
		
		chivertex = chivertex - nv;
		delete chivertex;
		
		return(0);
	}
	
	int InitialiseQuantities()
	{
		readbin ("Sigma", Sigmaloc-nvin, 2*nvin);
		
		calcEk(nk, Ek);
		
		calcGreal( Sigmacor, Sigmaloc, Ek, nv, nk, mu, beta, Gk);
		
		return(0);
	}
		
	int CalcChiBubble()
	{
		
		calcBubbleChi(Gk, beta, nv, nk, qtok, ksum, chibubble);
		return(0);
	}
	
	int CalcChiVertex(int FlocSep, int nvFloc)
	{
		calcVertexChi(Gk, Fup, Fdown, Flocup, Flocdown, FlocSep, beta, nv, nk, nvFloc, qtok, ksum, chivertex);
		return(0);
	}
	
	int WriteChis()
	{
		writebin("ChiBubble", chibubble-nv, (2*nv+1)*ndistk);
		writebin("ChiVertex", chivertex-nv, (2*nv+1)*ndistk);
		return(0);
	}
	
	int WriteGreal()
	{
		writebin("Greal", Gk - nv, 2*nv*nk*nk);
		return(0);
	}
	
	
};
