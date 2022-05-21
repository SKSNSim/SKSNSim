// File: FluxCalculation.cc
// Created  at 10 DEC 2021 by S.IZUMIYAMA
// Modified at 21 APR 2022 by S.IZUMIYAMA toward merge into SKSNSim
//
// Description:
// Class to load flux data file having data record like csv format
// -> #Energy flux

#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>
#include "FluxCalculation.hh"

#ifdef RETURNEXCEPATIONS
#include <stdexcept>
#endif

FluxCalculation::FluxCalculation()
{
      ene_flux_v = new std::vector<std::pair<double,double> >();
      lower_energy_bin_width = 0.1;
      upper_energy_bin_width = 0.1;
      return;
}

FluxCalculation::FluxCalculation(std::string fname)
{
      ene_flux_v = new std::vector<std::pair<double,double> >();
      this->loadFile(fname);
      return;
}

FluxCalculation::~FluxCalculation()
{
  delete ene_flux_v;
}

void FluxCalculation::loadFile(const std::string fname)
{
  ene_flux_v->clear();
  std::ifstream datafile(fname.c_str());
  double ene, flux;
  for(std::string line; std::getline(datafile, line); ){
    std::stringstream ss(line);
    ss >> ene >> flux;
    ene_flux_v->push_back(std::make_pair(ene, flux));
  }
  datafile.close();
  sortByEnergy();
  return;
}

void FluxCalculation::dumpFlux(std::ostream &out) const 
{
  for(std::vector<std::pair<double, double> >::iterator it = ene_flux_v->begin(); it != ene_flux_v->end(); it++)
    out << it->first << " MeV, " << it->second << std::endl;
}

void FluxCalculation::sortByEnergy()
{
  if(ene_flux_v == NULL) return;
  // for std::pair, first comparing first-component then comparing second-component
  std::sort(ene_flux_v->begin(), ene_flux_v->end());
  if(ene_flux_v->size() > 1){
    lower_energy_bin_width = (ene_flux_v->begin() + 1)->first - ene_flux_v->begin()->first;
    upper_energy_bin_width = (ene_flux_v->end()-1)->first - (ene_flux_v->end()-2)->first;
  }
  return;
}

double FluxCalculation::getFlux(const double nu_ene_MeV) const 
{
  // Assumed the data field ene_flux_v is sorted as lowest energy on first 
  // Calculate flux with linear interpolation
  const static double ERROR_CODE = -9999.;
  const static double OUTOFRANGE = -9998.;
  if(ene_flux_v == NULL) return ERROR_CODE;
  std::pair<double,double> bin = std::make_pair(-1, 0);
  std::pair<double,double> nextbin;
  for(std::vector<std::pair<double,double> >::iterator it = ene_flux_v->begin(); it+1 != ene_flux_v->end(); it++){
    // IMPROVEMENT: if you implement this with binary tree search, it will be faster
    if(it->first <= nu_ene_MeV && nu_ene_MeV < (it+1)->first){
      bin = *it;
      nextbin = *(it+1);
      break;
    }
  }
  if(bin.first == -1) {
#ifdef RETURNEXCEPATIONS
    throw  std::out_of_range("energy is out of range");
#endif
    return OUTOFRANGE;
  }

  double nuFlux = (nextbin.second - bin.second) * (nu_ene_MeV - bin.first) / (nextbin.first - bin.first) + bin.second;
#ifdef DEBUG
  std::cout << "FluxCalculation: nuEne = " << nu_ene_MeV <<", flux = " << nuFlux << ", bin = (" << bin.first << ", " << bin.second << ")" << std::endl;
#endif

  return nuFlux;
}

