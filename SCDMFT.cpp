using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <complex>
#include "routines.cpp"

// Self contained DMFT based on DF Sigma correction results

template <typename cpltype, typename fltype>
int PrecalciNu( cpltype* const iNu , const fltype beta, int nv)
//
// precalculates a vector of Matsubara frequencies for use in some functions.
{
	int i;
	fltype pibeta = acos( fltype(-1.))/beta;
	
	for(i=0;i< nv;i++)
	{
		*(iNu + i) = cpltype( 0. , pibeta * (2*i + 1) );
	}
	
	return (0);
}


int main(int argc, char* argv[])
{
	const double beta = atof(argv[1]); // Inverse Temperature
	const double U = atof(argv[2]);
	double mu = atof(argv[3]);
	const double t = 0.25;
	double Ef = atof(argv[4]);
	double p1 = atof(argv[5]);
	const int nk = atoi(argv[6]); // number of k-values per dimension
	const int nkin = atoi(argv[7]); // number of k-values per dimension in sigma-input
	const int nv = atoi(argv[8]); //Number of matsubara frequencies used
	const int nvin = atoi(argv[9]); //Number of matsubara frequencies in sigma-input
	const int niter = atoi(argv[10]); //Number of DMFT iterations (Usually ~100 is absolutely sufficient and cheap)
	
 	int i,j,k,l;
	int iin,jin;
	
	//Initialisation of DMFT-Variables:
	dcomp* const MatsuSigma = new dcomp [nv]; 
	dcomp* const MatsuDelta = new dcomp [nv]; 
	dcomp* const MatsuGloc = new dcomp [nv]; 
	dcomp* const MatsuG0 = new dcomp [nv]; 
	dcomp* const iNu = new dcomp [nv]; 
	
	dcomp* SigmaCorIn = new dcomp [2*nvin*nkin*nkin]; 
	dcomp* SigmaCor = new dcomp [nv*nk*nk]; 
	
	readbin ("SigmaCork" , SigmaCorIn , 2*nvin*nkin*nkin );
	
	SigmaCorIn += nvin; // set natural position of Sigma correction input.
	{
		const double coordtrans = double(nkin)/double(nk);
		double wi, wj;
		
		for(k=0; k<nvin; k++)
		{
			for(i = 0; i < nk; i++)
			{
				wi = coordtrans*i;
				iin = int(wi);
				wi -= iin;
				
				dcomp* duma;
				dcomp* dumb;
				
				duma = SigmaCorIn + k + iin*nkin*2*nvin;
				dumb = SigmaCorIn + k + ((iin+1)%nkin)*nkin*2*nvin;
				
				for(j = 0; j < nk; j++)
				{
					wj = coordtrans*j;
					jin = int(wj);
					wj -= jin;
					
					SigmaCor[k*nk*nk+i*nk+j] = (1.-wi)*((1.-wj)*duma[jin*2*nvin] + wj*duma[((jin+1)%nkin)*2*nvin] ) + wi*((1.-wj)*dumb[jin*2*nvin] + wj*dumb[((jin+1)%nkin)*2*nvin] );
				}
			}
		}
		for(k=nvin; k<nv; k++)
		{
			for(i = 0; i < nk; i++)
			{
				for(j = 0; j < nk; j++)
				{
					SigmaCor[k*nk*nk+i*nk+j] = 0.;
				}
			}
		}
		
	}
	cout << "reading of correction data done" << endl;
	
	
	double* Epsilons = new double [nk*nk]; // kinetic energies
	
	{
		double pi = acos( (-1.));
		double pistep = 2.*pi / nk;
		double epsa;
		
		for(i=0;i<nk;i++)
		{
			epsa = cos(i*pistep - pi);
			
			for(j=0;j<nk;j++)
			{
				*(Epsilons + i*nk + j) = -2. * t * (epsa + cos(j*pistep - pi));
			}
		}
	}
	
	i = PrecalciNu(iNu , beta , nv);
	
	for (i = 0; i<nv; i++) { MatsuSigma[i] = dcomp (p1*U, 0.);}  //Starting guess for Sigmaimag
	
	cout << "starting DMFT" << endl;
	
	{
		dcomp Gdum;
		dcomp dummy;
		dcomp dumsum;
		dcomp* dumpoint;
		double pU = p1 * U;
		double Up = U - pU;
		
		for (i = 0; i< niter; i++)
		{
			for (j = 0; j < nv; j++)
			{
				Gdum = 0.;
				dummy = *(iNu + j) + mu - *(MatsuSigma+j);
				dumpoint = SigmaCor + j*nk*nk;
				
				for (k = 0; k < nk; k++)
				{
					dumsum = 0.;
					for (l = 0; l < nk; l++)
					{
						dumsum += 1. / (dummy - *(Epsilons + k*nk + l) - *(dumpoint + k*nk + l) );
					}
					Gdum += dumsum;
				}
				*(MatsuGloc + j) = Gdum/double(nk*nk);
				*(MatsuG0 + j) = *(MatsuGloc + j) / (1. + *(MatsuGloc + j) * (*(MatsuSigma+j)) ) ;
				*(MatsuSigma+j) =  pU / (1. - Up * (*(MatsuG0 + j) ) ) ;
			}
		}
	}
	
	
	double occu = 0.;
	for (i = 0; i<nv; i++) 
	{ 
		occu += real(MatsuGloc[i]); 
	}
	
	cout << "DMFT done, writing Output" <<  endl;
	occu = 2. * occu / beta + 0.5;
	cout << "Occupation: " << occu << endl;
	
	ofstream outfile;
	//Output-scope
	{
	outfile.open("Gloc.txt" , ios::trunc );
	for (i = 0; i < nv; i++){
		outfile << imag(*(iNu + i)) << "   " << real(*(MatsuGloc + i)) << "     " << imag(*(MatsuGloc + i)) << "\n" ;
	}
	outfile.close();
	outfile.open("DMFTSigma.txt" , ios::trunc );
	for (i = 0; i < nv; i++){
		outfile << imag(*(iNu + i)) << "   " << real(*(MatsuSigma + i)) << "     " << imag(*(MatsuSigma + i)) << "\n" ;
	}
	outfile.close();
	
	outfile.open("Params.txt" , ios::trunc );
	outfile << "beta/ T = " << beta << " / " << 1. / beta << "\n" ;
	outfile << "U = " << U  << "\n" ;
	outfile << "mu = " << mu  << "\n" ;
	outfile << "t = " << t  << "\n" ;
	outfile << "Ef = " << Ef  << "\n" ;
	outfile << "p1 = " << p1  << "\n" ;
	outfile << "nk = " << nk  << "\n" ;
	outfile << "nv = " << nv  << "\n" ;
	outfile << "c-occupation = " << occu  << "\n" ;
	outfile.close();
	}
	
	
	dcomp* F = new dcomp[nv * nv * 4];
	F += nv * (2*nv+1);
	
	dcomp* Matsua = new dcomp [2 * nv]; 
	Matsua = Matsua + nv;
	const double ppUs = sqrt(p1 * (1.-p1) ) * U;
	
	for (i = 0; i < nv; i++){
		Matsua[i] = (MatsuSigma[i]) * (MatsuSigma[i] - U) / ppUs;
	}
	for (i = 0; i < nv; i++){
		Matsua[-i-1] = conj(Matsua[i]);
	}
	
	for (i = -nv; i < nv; i++){
            for (j = -nv; j < nv; j++){
                    F[i * 2*nv + j] = Matsua[i] * Matsua[j];
            }
	}
	Matsua = Matsua - nv;
	delete Matsua;
	
	
	outfile.open("FlocUp.txt" , ios::trunc );
	for (i = -nv; i < nv; i++){
            for (j = -nv; j < nv; j++){
                    outfile << 2 * (i) + 1 << "   " << 2*(j) + 1 << "   " << real(F[i * (2*nv) + j]) << "     " << imag(F[i * (2*nv) + j]) << "\n" ;
            }
	}
	outfile.close();
	
	
	F -= nv * (2*nv+1);
	
	writebin ("LambdaUp" , F , 4*nv*nv );
	writebin ("LambdaDown" , F , 4*nv*nv );
	
	delete F;
	
	dcomp* SuperG = new dcomp [2 * nv]; 
	SuperG = SuperG + nv;
	dcomp* SuperSigma = new dcomp [2 * nv]; 
	SuperSigma = SuperSigma + nv;
	
	for (i = 0; i < nv; i++){
		SuperG[i] = MatsuGloc[i];
		SuperSigma[i] = MatsuSigma[i];
	}
	for (i = 0; i < nv; i++){
		SuperG[-i-1] = conj(SuperG[i]);
		SuperSigma[-i-1] = conj(SuperSigma[i]);
	}
	
	//Prepare bare dual propagator
	
	double* inEpsilons = new double [nkin*nkin]; // kinetic energies
	
	{
		double pi = acos( (-1.));
		double pistep = 2.*pi / nkin;
		double epsa;
		
		for(i=0;i<nkin;i++)
		{
			epsa = cos(i*pistep - pi);
			
			for(j=0;j<nkin;j++)
			{
				*(inEpsilons + i*nkin + j) = -2. * t * (epsa + cos(j*pistep - pi));
			}
		}
	}
	
	dcomp* Gdual = new dcomp[2*nvin*nkin*nkin];
	Gdual += nvin;
	
	{
		dcomp* dumpoint;
		dcomp dummy;
		for (j = 0; j < nvin; j++)
		{
			dummy = *(iNu + j) + mu - *(MatsuSigma+j);
			dumpoint = SigmaCorIn + j;
			
			for (k = 0; k < nkin; k++)
			{
				for (l = 0; l < nkin; l++)
				{
// 					Gdual[k*nkin*2*nvin + l*2*nvin + j] = (1. / (dummy - *(inEpsilons + k*nkin + l) - *(dumpoint + (k*nkin + l) *2*nvin ) ) ) - SuperG[j];
					Gdual[k*nkin*2*nvin + l*2*nvin + j] = (1. / (dummy - *(inEpsilons + k*nkin + l) )) - SuperG[j];
					Gdual[k*nkin*2*nvin + l*2*nvin -j-1] = conj(Gdual[k*nkin*2*nvin + l*2*nvin + j]);
				}
			}
		}
	}
	
	Gdual -= nvin;
	writebin ("G0dual" , Gdual , 2*nvin*nkin*nkin );
	cout << nvin << " " << nkin*nkin << endl;
	delete Gdual;
	
	SuperG = SuperG - nv;
	SuperSigma = SuperSigma - nv;
	
	writebin ("G1" , SuperG , 2*nv );
	writebin ("Sigma" , SuperSigma , 2*nv );
	
	
	SigmaCorIn -= nvin; // reset natural position of Sigma correction input.
 	delete SigmaCorIn;

	delete SuperG;
	delete SuperSigma;
	delete MatsuSigma;
	delete MatsuDelta;
	delete MatsuGloc;
	delete MatsuG0;
	delete Epsilons;
	delete inEpsilons;
	return 0;
}
