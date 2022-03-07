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

const double nuElaEneMin = 0.0; //MeV, minimum neutrino energy 
const double nuElaEneMax = 150.0; //MeV, maximum neutrino energy 
const int nuElaEneNBins = 15000; // number of bins for nu energy
double nuElaEneBinSize = ( nuElaEneMax - nuElaEneMin ) / ( double )nuElaEneNBins;

class VectGenNuCrosssection
{
public:
  VectGenNuCrosssection();

  double CsNuebP_VB(double, double);
  double CsNuebP_SV(double, double);
  void DcsNuebP_VB(double, double, double&, double&);
  void DcsNuebP_SV(double, double, double&, double&);

  void ReadCsNuElastic();
  double CsNuElastic(int, double);

private:
  TFile* fCsElaFile;
  TTree* fCsElaTree;
  Double_t NuElaEnergy, CsElaNue, CsElaNux, CsElaNeb, CsElaNxb;

};

#endif
