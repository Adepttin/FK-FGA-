// Prepare Dispersion relation operator for nk*nk k-mesh
// Nearest neighbour-hopping square lattice is hardcoded so far
template <typename numbertype>
int calcEk(const int nk, numbertype* const Ek)
{
	int i,j;
	numbertype Ex;
	numbertype pi = acos( (-1.));
	numbertype pistep = 2.*pi / nk;
	const double t = 0.25;
	
	for(i = 0; i < nk; i++)
	{
		Ex = - 2. * t * cos( i * pistep - pi );
		for(j = 0; j < nk; j++)
		{
			*(Ek + i*nk + j) = Ex - 2. * t * cos( j * pistep - pi );
		}
	}
	
	return(0);
}


// Prepare velocity-x-operator for nk*nk k-mesh
// Nearest neighbour-hopping square lattice is hardcoded so far
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

// Prepare fermionic Matsubara frequencies
template <typename cpltype, typename fltype>
int PrecalciNu( cpltype* const iNu , const fltype beta, int nv)
//
// precalculates a vector of Matsubara frequencies for use in some functions.
{
	int i;
	fltype pibeta = acos( fltype(-1.))/beta;
	
	for(i=-nv;i< nv;i++)
	{
		*(iNu + i) = cpltype( 0. , pibeta * (2*i + 1) );
	}
	
	return (0);
}

// calculate Bubble-contribution to conductivity from Kubo's formula
// Use Gk (already dressed by non-local Sigma-corrections), Matsubara frequencies and x-velocities vx as input
// Correct for finite frequency box by subtracting the trivial propagators 1/(iv * (iv + iw)) which give 0.25 for w = 0 and 0 otherwise when summed over frequencies
// write result to ohm (conductivity times omega!)
template <typename numbertype , typename fltype>
int calcohm(const numbertype* const Gk, const numbertype* const Matsus, const fltype* const vx, const int nv2, const int nw2, const int nk, const fltype beta, numbertype* const ohm )
{
	int i,j,k,l;
	numbertype dummy, dumg;
	const int norm = nk*nk;
	numbertype duma, dumb;
	const fltype flnorm = norm;
	
	fltype fldummy, zeroohm;
	
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
	
	
	for(k = 0; k < nw2+1; k++)
	{
		*(ohm + k) = 0.;
		
		//dummy is a primitive anti-absorbtion measure for numerical summation
		for(j = 0; j < nk; j++)
		{
			
			dummy = 0.;
			for(l = 0; l < nk; l++)
			{
				for(i = -nv2; i < nv2-k; i++)
				{
					duma = *(Gk + (j*nk + l)*2*nv2 + i);
					dumb = *(Gk + (j*nk + l)*2*nv2 + i+k);
					dumg = 1./(Matsus[i] * Matsus[i+k]); 
					dummy += ( duma * dumb - dumg ) * *(vx + j*nk + l) * *(vx + j*nk + l);
				}
				
			}
			*(ohm + k) += dummy;
		}
		
		*(ohm + k) = *(ohm + k) / (flnorm);
		*(ohm + k) = *(ohm + k) / (beta*beta);
	}
	
	*(ohm) += zeroohm;
	
	for(k = 1; k < nw2+1; k++)
	{
		*(ohm - k) = *(ohm + k);
	}
	
	
	return(0);
}


// calculate Vertex-contribution to conductivity from Kubo's formula
// Use Fq, Gk (already dressed by non-local Sigma-corrections), Matsubara frequencies and x-velocities vx as input
// Correct for finite frequency box by subtracting the trivial propagators C * 1/(iv * (iv + iw))^2 which give w-constant contributions per kk' set
// write result to ohm (conductivity times omega!)
template <typename numbertype , typename fltype>
int calcconohm(const numbertype* const* const Fqdown, const numbertype* const Gk, const numbertype* const Matsus, const fltype* const vx, const int nv2, const int nw2, const int nk, const fltype beta, numbertype* const ohm )
{
	int i,k,kp,v;
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
	
	for(i = 0; i < nw2 + 1; i++)
	{
		*(ohm + i) = 0.;
		for(k = 0; k < norm; k++)
		{
			dummya = 0.;
			
			for(kp = 0; kp < norm; kp++)
			{
				dummyb = 0.;
				
				Casym = Fqdown[k*norm + kp][(nv2-1)*(2*nv2)-nv2];
				
				for(v = -nv2; v < nv2-i; v++)
				{
					dumg = 1./(Matsus[v] * Matsus[i+v]); 
					dumg = dumg*dumg;
					dumg2 = Gk[k*2*nv2 + v] * Gk[kp*2*nv2 + v] * Gk[kp*2*nv2 + (v+i)] * Gk[k*2*nv2 + (v+i)];
					
					dummyb += ( (Fqdown[k*norm + kp][(v)*(2*nv2+1)+i] * dumg2) - (Casym * dumg) ) * vx[k] * vx[kp];
				}
				
				dummya += (dummyb/(beta*beta) + Casym * wref[i] * vx[k] * vx[kp]);
			}
			
			*(ohm + i) += dummya/(flnorm*flnorm);
		}
	}
	
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
	int nk, nv, nvin;
	double beta,mu;
	double* vx;
	double* Ek;
	
	//one-particle quantities
	dcomp* Sigmaloc;
	dcomp* Sigmacor;
	dcomp* Matsus;
	dcomp* Gk;
	
	//conductivities
	dcomp* ohmbubble;
	dcomp* ohmvertex;
	
	//Vertex quantities
	dcomp ** Fdown;
	
	template <class Parquet>
	ConductivityObject(Parquet MyParquet, double inbeta, double inmu)
	{
		nk = MyParquet.nk;
		nv = MyParquet.nv;
		nvin = MyParquet.nvin;
		
		Sigmacor = MyParquet.Sigmacor;
		Fdown = (MyParquet.Fdown)[MyParquet.ndistk - 1];
		
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
// 		cout << "Greal calculated" << endl;
		return(0);
	}
	
	int CalcCondBubble()
	{
		calcohm(Gk, Matsus, vx, nv, nv, nk, beta, ohmbubble);
		return(0);
	}
	
	int CalcCondVertex()
	{
		calcconohm(Fdown, Gk, Matsus, vx, nv, nv, nk, beta, ohmvertex );
		cout << "Vertex calculated" << endl;
		return(0);
	}
	
	int WriteConductivities()
	{
		writebin("BubbleOhm", ohmbubble-nv, 2*nv+1);
		writebin("VertexOhm", ohmvertex-nv, 2*nv+1);
		
		return(0);
	}
	
};

