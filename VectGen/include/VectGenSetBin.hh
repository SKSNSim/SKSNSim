#ifndef VECTGENSETBIN_H_INCLUDED
#define VECTGENSETBIN_H_INCLUDED

//double eEneThr = 3.0; // MeV, electron energy threshold

//const double nuEneMin = 3.0; //MeV, minimum neutrino energy 
//const double nuEneMax = 80.0; //MeV, maximum neutrino energy 
//const int nuEneNBins = 770; // number of bins for nu energy
//double eEneThr = 1.0;// MeV, electron total energy threshold
double eEneThr = 5.0;// MeV, electron total energy threshold

const double nuEneMin = 0.0; //MeV, minimum neutrino energy 
const double nuEneMax = 100.0; //MeV, maximum neutrino energy 
const int nuEneNBins = 1000; // number of bins for nu energy
double nuEneBinSize;

const double tStart = 0.;    //sec, SN explosion start time
const double tEnd = 20.0;    //sec, SN explotion end time
const double tBinSize = 1.0e-3;  //sec, time interval
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
