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
// Self contained DF based on DMFT results


int main(int argc, char* argv[])
{
	const double beta = atof(argv[1]); // Inverse Temperature
	const double U = atof(argv[2]);
	double mu = atof(argv[3]);
	const double t = 0.25;
	double Ef = atof(argv[4]);
	double p1 = atof(argv[5]);
	const int nk = atoi(argv[6]); // number of k-values per dimension
	const int nkin = atoi(argv[7]); // number of k-values per dimension in input
	const int nv = atoi(argv[8]); //Number of matsubara frequencies used
	const int nvin = atoi(argv[9]); //Number of matsubara frequencies in input
	const double relconv = atoi(argv[10]); //desired relative convergence for DF-corrections
	const int maxit = atoi(argv[11]); //maximum number of iterative DF loops
	
	int i;
	
	cout << "Blubb" <<  endl;
	DFParquetParams ParqObj = DFParquetParams (nk , nv , nkin , nvin);
	
	cout << "Object built" <<  endl;
	
	cout << "Initialising 1P" <<  endl;
	ParqObj.InitialiseOneParticle();
	cout << "Initialising Kuantities" <<  endl;
	ParqObj.InitialiseKQuantities();
	cout << "Initialising V" <<  endl;
	ParqObj.InitialiseVertexStorage();
	
	cout << "Reading DMFT" <<  endl;
	ParqObj.ReadDMFT();
	
	//reset phi and initialise F and Gamma to Floc
	ParqObj.ResetVertex();
	ParqObj.Parquetiter();
	
	
	//reading Sigmadual and updating Gdual
	ParqObj.ReadDualSig();
	ParqObj.UpdateGdual();
	
	for(i = 0; i < maxit; i++)
	{
		cout << "BS Iteration " << i << "/" << maxit << endl;
		ParqObj.BSiter();
		cout << "Parquet" << endl;
		ParqObj.Parquetiter();
		
		cout << "Calculating Sdual" <<  endl;
		ParqObj.SigCalc();
		ParqObj.UpdateGdual();
	}
	

	
	ConductivityObject CondObj = ConductivityObject(ParqObj, beta, mu);
	CondObj.InitialiseStorage();
	CondObj.InitialiseQuantities();
	
	CondObj.CalcCondBubble();
	CondObj.CalcCondVertex();
	
	CondObj.WriteConductivities();
	
	CondObj.DeleteStorage();
	
	cout << "Writing Sdual" <<  endl;
	ParqObj.FlexDualToRealSig(0);
	ParqObj.WriteSigCors();
	ParqObj.WriteDualSig();
	ParqObj.WriteGdual();
	
	cout << "Deleting Khelper" <<  endl;
	ParqObj.DeleteKQuantities();
	cout << "Writing 1P quantities" <<  endl;
	ParqObj.DeleteOneParticle();
	cout << "Deleting Vertices" <<  endl;
	ParqObj.DeleteVertexStorage();
	
	ofstream outfile;
	//Output-scope
	{
	outfile.open("DFParquetParams.txt" , ios::trunc );
	outfile << "beta/ T = " << beta << " / " << 1. / beta << "\n" ;
	outfile << "U = " << U  << "\n" ;
	outfile << "mu = " << mu  << "\n" ;
	outfile << "t = " << t  << "\n" ;
	outfile << "Ef = " << Ef  << "\n" ;
	outfile << "p1 = " << p1  << "\n" ;
	outfile << "nk = " << nk  << "\n" ;
	outfile << "nv = " << nv  << "\n" ;
	outfile.close();
	}
	
	
	return 0;
}
