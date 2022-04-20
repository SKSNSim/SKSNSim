// File: FluxCalculation.hh
// Created  at 10 DEC 2021 by S.IZUMIYAMA
// Modified at 21 APR 2022 by S.IZUMIYAMA towards merge into SKSNSim

#ifndef _FLUXCALCULATION_HH_INCLUDED_
#define _FLUXCALCULATION_HH_INCLUDED_

#include <string>
#include <vector>
#include <iostream>

// #define RETURNEXCEPATIONS 
// This enables throwing exceptions to notify some errror (e.g. the specified energy is out of range)
#ifdef RETURNEXCEPATIONS
#include <typeinfo>
#endif

class FluxCalculation {
  private:
    std::vector<std::pair<double,double> > *ene_flux_v = NULL; /* size_t -> <energy, flux> */
    double lower_energy_bin_width;
    double upper_energy_bin_width;
    void sortByEnergy();

  public:
    const static size_t NBIN = 500;
    FluxCalculation();
    FluxCalculation(const std::string /*fname*/);
    ~FluxCalculation();
    void loadFile(const std::string /*fname*/);
    double getFlux(const double /* nu energy in MeV */) const;
    inline double getFluxLimit(const bool lower_limit = true /* if false -> return upper limit*/) const{
#ifdef RETURNEXCEPATIONS
      if ( ene_flux_v == NULL ) throw std::bad_typeid("ene_flux_v is not allocated");
#endif
      if ( ene_flux_v == NULL ) return -9999;
      return  lower_limit? ene_flux_v->front().first - lower_energy_bin_width/2.0: ene_flux_v->back().first + upper_energy_bin_width/2.0;
    }
    // This throws std::out_of_range when the energy is out of range of loaded flux data
    void dumpFlux(std::ostream &out = std::cout) const;
    inline double getBinnedFlux(const int bin) const {return ene_flux_v->at(bin).second;};
    inline double getBinnedEnergy(const int bin) const {return ene_flux_v->at(bin).first;};
    inline size_t getNBins() const {return (ene_flux_v != NULL) ? ene_flux_v->size(): 0;};
};

#endif
