using namespace std;
// #include <stdio.h>
// #include <stdlib.h>
// #include <fstream>
// #include <iostream>
// #include <assert.h>
// #include <cmath>
// #include <complex>

// #include "routines.cpp"
#include "DFphLadderObj.cpp"
// #include "DFConductivity.cpp"
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
	
//  	int i,j,k,l;
	
	cout << "Blubb" <<  endl;
	DFphParams LadderObj = DFphParams (nk , nv , nkin , nvin);
	
	cout << "Object built" <<  endl;
	
	cout << "Initialising 1P" <<  endl;
	LadderObj.InitialiseOneParticle();
	cout << "Initialising Kuantities" <<  endl;
	LadderObj.InitialiseKQuantities();
	cout << "Initialising V" <<  endl;
	LadderObj.InitialiseVertexStorage();
	
	cout << "Reading DMFT" <<  endl;
	LadderObj.ReadDMFT();

//reading Sigmadual and updating Gdual
	LadderObj.ReadDualSig();
	LadderObj.UpdateGdual();
	
	for(i = 0; i < maxit; i++)
	{
		//BS only for ph
		cout << "BS Iteration " << i << "/" << maxit << endl;
		LadderObj.LadderCalc();
		//Parquet and calculation of Sigma
		cout << "Sigma calculation" << endl;
		LadderObj.SigCalc();
		LadderObj.UpdateGdual();
	}
	
	cout << "Update Gdual" <<  endl;
	//ParqObj.SigCalc();
	LadderObj.UpdateGdual();
	
	cout << "Writing Sdual" <<  endl;
	LadderObj.FlexDualToRealSig(0);
	LadderObj.WriteSigCors();
	LadderObj.WriteDualSig();
	LadderObj.WriteGdual();
	
    ConductivityObject CondObj = ConductivityObject(LadderObj, beta, mu);
	cout << "Initialise Conductivity storage" <<  endl;
	CondObj.InitialiseStorage();
	CondObj.InitialiseQuantities();
	
	cout << "Calculate Bubble" <<  endl;
	CondObj.CalcCondBubble();
	cout << "Calculate ConnOhm" <<  endl;
	CondObj.CalcCondVertex();
	
	CondObj.WriteConductivities();
	
	CondObj.DeleteStorage();
// 	ConductivityObject CondObj = ConductivityObject(ParqObj, beta, mu);
// 	CondObj.InitialiseStorage();
// 	CondObj.InitialiseQuantities();
// 	
// 	CondObj.CalcCondBubble();
// 	CondObj.CalcCondVertex();
// 	
// 	CondObj.WriteConductivities();
// 	
// 	CondObj.DeleteStorage();
	
	
	cout << "Deleting Khelper" <<  endl;
	LadderObj.DeleteKQuantities();
	cout << "Writing 1P quantities" <<  endl;
	LadderObj.DeleteOneParticle();
	cout << "Deleting Vertices" <<  endl;
	LadderObj.DeleteVertexStorage();
	
	ofstream outfile;
	//Output-scope
	{
	outfile.open("DFphLadderParams.txt" , ios::trunc );
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
