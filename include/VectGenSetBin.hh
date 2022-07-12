#ifndef VECTGENSETBIN_H_INCLUDED
#define VECTGENSETBIN_H_INCLUDED

//double eEneThr = 3.0; // MeV, electron energy threshold

//const double nuEneMin = 3.0; //MeV, minimum neutrino energy 
//const double nuEneMax = 80.0; //MeV, maximum neutrino energy 
//const int nuEneNBins = 770; // number of bins for nu energy
//double eEneThr = 0.52;// MeV, electron total energy threshold
//double eEneThr = 1.0;// MeV, electron total energy threshold
double eEneThr = 5.0;// MeV, electron total energy threshold

//const double nuEneMin = 0.0; //MeV, minimum neutrino energy 
//const double nuEneMax = 100.0; //MeV, maximum neutrino energy 
//const double nuEneMax = 60.0; //Nakahata
const double nuEneMin = 0.0; //MeV, minimum neutrino energy (Horiuchi SRN)
const double nuEneMax = 300.0; //MeV, maximum neutrino energy (Horiuchi SRN)
const int nuEneNBins = 3000; // number of bins for nu energy
//const int nuEneNBins = 120; // Nakahata
double nuEneBinSize;

const double tStart = 0.;    //sec, SN explosion start time
//const double tStart = 0.02;    //Nakahata
const double tEnd = 20.0;    //sec, SN explotion end time
//const double tEnd = 18.0;    //sec, SN explotion end time
const double tBinSize = 1.0e-3;  //sec, time interval
//const double tBinSize = 1.0e-2;  //Nakahata
int tNBins;

const double costMin = -1.;
const double costMax = 1.;
const int costNBins = 1000;
double costBinSize;


void VectGenSetBinValues( void )
{
	// set bin
	nuEneBinSize = ( nuEneMax - nuEneMin ) / ( double )nuEneNBins;

	tNBins = ( int ) ( ( tEnd - tStart ) / tBinSize );
	//tNBins -= 2; //to prevent nan

	costBinSize = ( costMax - costMin ) / ( double )costNBins;
}

#endif // SETBIN_H_INCLUDED
