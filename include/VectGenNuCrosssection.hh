/**
 * @file VectGenNuCrosssection.hh
 *
 * @date 2017-12-06
 * @author Y.Koshio
 */

#ifndef VECTGENNUCROSSSECTION_HH
#define VECTGENNUCROSSSECTION_HH
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

/**
 * @class VectGenNuCrosssection
 */

extern "C" {
	double sl_nue_dif_rad_( double *, double * );
	double sl_neb_dif_rad_( double *, double * );
	double sl_num_dif_rad_( double *, double * );
	double sl_nmb_dif_rad_( double *, double * );
}

const double nuElaEneMin = 0.0; //MeV, minimum neutrino energy 
const double nuElaEneMax = 150.0; //MeV, maximum neutrino energy 
const int nuElaEneNBins = 15000; // number of bins for nu energy
double nuElaEneBinSize = ( nuElaEneMax - nuElaEneMin ) / ( double )nuElaEneNBins;

class VectGenNuCrosssection
{
public:
  VectGenNuCrosssection();

  double CsNuebP_VB(double);
  double CsNuebP_SV(double);
  void DcsNuebP_VB(double, double, double&, double&);
  void DcsNuebP_SV(double, double, double&, double&);

  void ReadCsNuElastic();
  double CsNuElastic(int, double, int);
  double calcElectronTotEnergyElastic(const double, const double);
  double calcDeEneDcostElastic(const double, const double);
  double calcCosTthElastic(const double, const double);
  double getEnu(double, double);

private:
  TFile* fCsElaFile;
  TTree* fCsElaTree;
  Double_t NuElaEnergy, CsElaNue, CsElaNux, CsElaNeb, CsElaNxb;

};

#endif
