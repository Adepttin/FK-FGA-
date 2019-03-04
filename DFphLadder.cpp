using namespace std;
//#include <stdio.h>
//#include <stdlib.h>
//#include <fstream>
//#include <iostream>
//#include <assert.h>
//#include <cmath>
//#include <complex>

//#include "routines.cpp"
#include "DFphLadderObj.cpp"

//DF ladder calculation based on DMFT results

int main(int argc, char* argv[])
{
	//read in parameters
	const double beta = atof(argv[1]); //inverse temperature
	const double U = atof(argv[2]); //interaction strength
	double mu = atof(argv[3]); //chemical potential
	const double t = 0.25; //hopping parameter
	double p1 = atof(argv[4]); //occupation of localized f-electrons
	const int nk = atoi(argv[5]); //number of k-values per dimension
	const int nkin = atoi(argv[6]); //number of k-values per dimension in DMFT input
	const int nv = atoi(argv[7]); //number of Matsubara frequencies used
	const int nvin = atoi(argv[8]); //number of Matsubara frequencies in DMFT input
	
	int i;
	
	//create instance of class DFppParams
	DFphParams LadderObj = DFphParams (nk, nv, nkin, nvin);
	cout << "Object built" <<  endl;
	
	//initialise quantities and allocate memory
	cout << "Initialising 1P" <<  endl;
	LadderObj.InitialiseOneParticle();
	cout << "Initialising Kuantities" <<  endl;
	LadderObj.InitialiseKQuantities();
	cout << "Initialising V" <<  endl;
	LadderObj.InitialiseVertexStorage();
	
	//read in DMFT results
	cout << "Reading DMFT" <<  endl;
	LadderObj.ReadDMFT();

	//read in dual self energy and update dual propagator
	LadderObj.ReadDualSig();
	LadderObj.UpdateGdual();
	
	//calculate ph ladder
	cout << "Calculate ph ladder" << endl;
	LadderObj.LadderCalc();
	
	//update dual self energy and propagator
	cout << "Calculate Sigma" << endl;
	LadderObj.SigCalc();
	LadderObj.UpdateGdual();
	
	//write self energy corrections, dual self energy and propagator
	cout << "Writing Sdual" <<  endl;
	LadderObj.FlexDualToRealSig(0);
	LadderObj.WriteSigCors();
	LadderObj.WriteDualSig();
	LadderObj.WriteGdual();
	
	//calculate current-current correlation function
	
	//create instance of class ConductivityObject	
	ConductivityObject CondObj = ConductivityObject(LadderObj, beta, mu);
	
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
	
	//release memory
	cout << "Deleting Khelper" <<  endl;
	LadderObj.DeleteKQuantities();
	cout << "Writing 1P quantities" <<  endl;
	LadderObj.DeleteOneParticle();
	cout << "Deleting Vertices" <<  endl;
	LadderObj.DeleteVertexStorage();
	
	//write parameters to txt file
	ofstream outfile;
	//Output-scope
	{
	outfile.open("DFphLadderParams.txt" , ios::trunc );
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
