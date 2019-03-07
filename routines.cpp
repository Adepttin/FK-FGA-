// list of routines which may be needed for more than one step of the 2-particle DF parquet

typedef complex<double> dcomp; //double precision complex type
typedef double dflt;           //double precision float type
const dcomp compli = dcomp(0.,1.); //complex i
const double pi = 4. * atan(1.); //pi
const double hopt = -0.25; //hopping parameter t

//reading binary data
template <typename numbertype>
int readbin (const char* filename , numbertype* const target , const int size )
{
	const int length = sizeof(numbertype)*size;
	ifstream in;
	cout << "reading " << length << " of data from " << filename << endl;
	
	in.open ( filename , ios::binary );
	if(in.is_open())
	{
		in.read ((char*)target,length);
		in.close();
	}
	else
	{
		cout << "unable to read from " << filename << endl;
		exit(EXIT_FAILURE);
	}
	
	return(0);
}

//writing binary data
template <typename numbertype>
int writebin (const char* filename , numbertype* const source , const int size )
{
	const int length = sizeof(numbertype)*size;
	ofstream out;
	cout << "writing " << length << " of data to " << filename << endl;
	
	out.open ( filename , ios::binary );
	if(out.is_open())
	{
		out.write ((char*)source,length);
		out.close();
	}
	else
	{
		cout << "unable to write to " << filename << endl;
		exit(EXIT_FAILURE);
	}
	
	return(0);
}

//prepare dispersion relation operator Ek for nk*nk k-mesh
//nearest neighbour-hopping square lattice is hardcoded so far
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

//calculate Green's function of real fermions, based on DF self energy corrections Sigmacor, DMFT self energy Sigma, dispersion relation epsilon, mu and beta 
//lattice size nk and number of Matsubara frequencies as in DF calculation
//result is written to Gk[v][kx][ky]
template <typename numbertype , typename fltype >
int calcGreal( const numbertype* const Sigmacor, const numbertype* const Sigma, const fltype* const epsilon, const int nv, const int nk, const fltype mu, const fltype beta, numbertype* const Gk)
{
	int i,j,k;
	numbertype dummy;
	numbertype* dummya;
	const numbertype* dummyb;
	
	for(i = 0; i < nk; i++)
	{
		for(j = 0; j < nk; j++)
		{
			dummya = Gk + (nk*i+j)*2*nv;
			dummyb = Sigmacor + (nk*i+j)*2*nv;
			for(k = -nv; k < nv; k++)
			{
				dummy = compli*((2 * k + 1)*pi/beta) + mu - Sigma[k];
				
				
				dummya[k] = 1./(dummy - dummyb[k] - epsilon[i*nk + j]);
			}
		}
	}
	
	return(0);
}

// // calculate q-shifted sum over quadratic (2D) Brillouin zone of Gk and Gkp = sum (Gk(k1) * Gkp(k1 + q)) 
// // calculation is carried out for (2*nv2 + boson) frequencies. Gk and Gkp have to be arranged as Gk[v1][kx][ky]
// // Result is written to chi[v]
// // Typically used to calculate bare susceptibility (w,q) with Gkp = Gk - w*nk*nk
// template <typename numbertype>
// int calchi(numbertype* const Gk, numbertype* const Gkp, const int nv2, const int nk, const int q1, const int q2, numbertype* const chi , const int boson = 0)
// {
// 	int i,j,k;
// 	numbertype dummy;
// 	const int effq1 = q1 + nk/2;
// 	const int effq2 = q2 + nk/2;
// 	const int norm = nk*nk;
// 	numbertype *duma, *dumb;
// 	
// 	for(k = -nv2; k < nv2+boson; k++)
// 	{
// 		*(chi + k) = 0.;
// 		for(i = 0; i < nk; i++)
// 		{
// 			duma = Gk + norm*k + i*nk;
// 			dumb = Gkp + norm*k + ((i+effq1)%nk)*nk;
// 			
// 			dummy = 0.;
// 			//dummy is a primitive anti-absorbtion measure for numerical summation
// 			for(j = 0; j < nk; j++)
// 			{
// 				dummy += *(duma + j) * *(dumb + ( (j+effq2)%nk ) );
// 				
// 			}
// 			
// 			*(chi + k) += dummy;
// 		}
// 		*(chi + k) = *(chi + k) / ((dflt)norm);
// 		
// 	}
// 	
// 	return(0);
// }
// 
// // calculate dual dispersion relation for quadratic nearest neigbour hopping with t=-0.25 
// // Result is written to epsilon[kx][ky]
// template <typename fltype >
// int calcNNepsilon(fltype* const epsilon, const int nk)
// {
// 	int i,j;
// 	const int nk2 = nk/2;
// 	fltype* coseps = new fltype[nk];                            //1D epsilons for k-mesh
// 	const fltype kstep = pi/nk2;
// 	fltype dumfloat;
// 	
// 	for(i = 0; i < nk; i++)
// 	{
// 		*(coseps+i) = 2 * hopt * cos( (i-nk2)*kstep );
// 	}
// 	
// 	for(i = 0; i < nk; i++)
// 	{
// 		dumfloat = coseps[i];
// 		
// 		for(j = 0; j < nk; j++)
// 		{
// 			epsilon[i*nk + j] = dumfloat + coseps[j];
// 		}
// 	}
// 	
// 	delete(coseps);
// 	
// 	return(0);
// }
// 
// // calculate dual Green's function, based on available dispersion relation, DMFT-sigma, local DMFT Green's function, mu and beta 
// // Result is written to Gk[v][kx][ky]
// template <typename numbertype , typename fltype >
// int calcGdual( const numbertype* const Sigma, const numbertype* const G1, const fltype* const epsilon, const int nv1, const int nk, const fltype mu, const fltype beta, numbertype* const Gk)
// {
// 	int i,j,k;
// 	numbertype dummy, ymmud;
// 	numbertype* dummya;
// 	
// 	for(k = -nv1; k < nv1; k++)
// 	{
// 		dummy = compli*((2 * k + 1)*pi/beta) + mu - Sigma[k];
// 		ymmud = G1[k];
// 		
// 		dummya = Gk + nk*nk*k;
// 		
// 		for(i = 0; i < nk; i++)
// 		{
// 			for(j = 0; j < nk; j++)
// 			{
// 				 dummya[i*nk + j] = 1./(dummy - epsilon[i*nk + j]) - ymmud;
// 			}
// 		}
// 	}
// 	
// 	return(0);
// }
// 
// 
// 
// 
// // dress dual Green's function, Gk[v][kx][ky] by dual Self energy from files 
// // Result is written to Gk[v][kx][ky]
// template <typename numbertype>
// int dressGdual( const int nv1, const int v2o , const int nk, numbertype* const Gk)
// {
// 	const int nv2o = v2o/2;
// 	int i,j,k;
// 	
// 	numbertype* dumsig = new numbertype [v2o];
// 	dumsig += nv2o;
// 	
// 	numbertype dummy, ymmud;
// 	numbertype* dummya;
// 	numbertype* dummyb;
// 	
// 	for(i = 0; i < nk; i++)
// 	{
// 		stringstream ss;
// 		
// 		for(j = 0; j <= i; j++)
// 		{
// 			ss.str(string());
// 			ss << "DFsig/Sig" <<"(" << i << ")(" << j << ")of(" << nk << ").out";
// 			readbin ( ss.str().c_str() , dumsig-nv2o , v2o );
// 			
// 			dummya = Gk + i*nk + j;
// 			dummyb = Gk + j*nk + i;
// 			
// 			for(k = -nv2o; k < nv2o; k++)
// 			{
// 				ymmud = *(dummya + k*nk*nk);
// 				dummy = ymmud/( 1. - ymmud * dumsig[k] );
// 				
// 				*(dummya + k*nk*nk) = dummy;
// 				*(dummyb + k*nk*nk) = dummy;
// 			}
// 		}
// 	}
// 	
// 	dumsig -= nv2o;
// 	delete dumsig;
// 	
// 	return(0);
// }
// 
// 
// // calculate dual Green's function, quadratic nearest neigbour hopping
// // Result is written to Gk[v][kx][ky]
// template <typename numbertype , typename fltype >
// int calcNNGdual( const numbertype* const Sigma, const numbertype* const G1, const int nv1, const int nk, const fltype mu, const fltype beta, numbertype* const Gk)
// {
// 	int retval = 0;
// 	fltype* epsilon = new fltype[nk*nk];
// 	retval += calcNNepsilon(epsilon, nk);
// 	retval += 2*calcGdual(Sigma, G1, epsilon, nv1, nk, mu, beta, Gk);
// 	
// 	delete(epsilon);
// 	
// 	return(retval);
// }



