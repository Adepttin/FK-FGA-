using namespace std;
// #include <stdio.h>
// #include <stdlib.h>
// #include <fstream>
// #include <iostream>
// #include <assert.h>
// #include <cmath>
// #include <complex>

// #include "routines.cpp"
#include "DFParquetObj.cpp"
#include "DFConductivity.cpp"
#include "DFMagnetism.cpp"

// Self contained DF parquet based on DMFT results

int main(int argc, char* argv[])
{
	//read in parameters
	const double beta = atof(argv[1]); //inverse Temperature
	const double U = atof(argv[2]); //interaction strength
	double mu = atof(argv[3]); //chemical potential
	const double t = 0.25; //hopping parameter
	double p1 = atof(argv[4]); //occupation of localized f-electrons
	const int nk = atoi(argv[5]); //number of k-values per dimension
	const int nkin = atoi(argv[6]); //number of k-values per dimension in DMFT input
	const int nv = atoi(argv[7]); //number of Matsubara frequencies used
	const int nvin = atoi(argv[8]); //number of Matsubara frequencies in DMFT input
	const int maxit = atoi(argv[9]); //number of iterative DF loops
	
	int i;
	
	//create instance of class DFParquetParams
	DFParquetParams ParqObj = DFParquetParams (nk , nv , nkin , nvin);
	cout << "Object built" <<  endl;
	
	//initialise quantities and allocate memory
	cout << "Initialising 1P" <<  endl;
	ParqObj.InitialiseOneParticle();
	cout << "Initialising Kuantities" <<  endl;
	ParqObj.InitialiseKQuantities();
	cout << "Initialising V" <<  endl;
	ParqObj.InitialiseVertexStorage();
	
	//read in DMFT results
	cout << "Reading DMFT" <<  endl;
	ParqObj.ReadDMFT();
	
	//reset vertices
	ParqObj.ResetVertex();	
	
	//read in dual self energy and update dual propagator
	ParqObj.ReadDualSig();
	ParqObj.UpdateGdual();
	
	//DF loop
	for(i = 0; i < maxit; i++)
	{
		//call Bethe-Salpeter equations
		cout << "BS Iteration " << i+1 << "/" << maxit << endl;
		ParqObj.BSiter();
		//call parquet equation
		cout << "Parquet" << endl;
		ParqObj.Parquetiter();
	}
	
	//update dual self energy and propagator
	cout << "Calculating Sdual" <<  endl;
	ParqObj.SigCalc();
	ParqObj.UpdateGdual();
	
	//write self energy corrections, dual self energy and propagator
	cout << "Writing Sdual" <<  endl;
	ParqObj.FlexDualToRealSig(0);
	ParqObj.WriteSigCors();
	ParqObj.WriteDualSig();
	ParqObj.WriteGdual();
	
	
	//calculate current-current correlation function
	
	//create instance of class ConductivityObject	
	ConductivityObject CondObj = ConductivityObject(ParqObj, beta, mu);
	
	//initialise quantities and allocate memory
	CondObj.InitialiseStorage();
	CondObj.InitialiseQuantities();
	
	//calculate bubble
	CondObj.CalcCondBubble();
	
	//calculate vertex corrections
	CondObj.CalcCondVertex();
	
	//write bubble and vertex corrections
	CondObj.WriteConductivities();
	
	//release memory
	CondObj.DeleteStorage();
	
	
	//calculate density-density correlation function

	//create instance of class MagnetismObject
	MagnetismObject MagObj = MagnetismObject(ParqObj, beta, mu);
	
	//initialise quantities and allocate memory
	MagObj.InitialiseStorage();
	MagObj.InitialiseQuantities();
	
	//calculate bubble
	MagObj.CalcChiBubble();

	//calculate vertex correction
	MagObj.CalcChiVertex();
	
	//write bubble and vertex corrections
	MagObj.WriteChis();
	
	//release memory
	MagObj.DeleteStorage();
	
	//release memory
	cout << "Deleting Khelper" <<  endl;
	ParqObj.DeleteKQuantities();
	cout << "Writing 1P quantities" <<  endl;
	ParqObj.DeleteOneParticle();
	cout << "Deleting Vertices" <<  endl;
	ParqObj.DeleteVertexStorage();
	
	//write parameters to txt file
	ofstream outfile;
	//Output-scope
	{
	outfile.open("DFParquetParams.txt" , ios::trunc );
	outfile << "beta/ T = " << beta << " / " << 1. / beta << "\n" ;
	outfile << "U = " << U  << "\n" ;
	outfile << "mu = " << mu  << "\n" ;
	outfile << "t = " << t  << "\n" ;
	outfile << "p1 = " << p1  << "\n" ;
	outfile << "nk = " << nk  << "\n" ;
	outfile << "nv = " << nv  << "\n" ;
	outfile.close();
	}
	
	
	return 0;
}
