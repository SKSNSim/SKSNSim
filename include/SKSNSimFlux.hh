/****************************************
 * File: SKSNSimFlux.hh
 * Description:
 *   Flux model for supernovae simulation
 ****************************************/


#ifndef __SKSNSIMFLUX_H_INCLUDED__
#define __SKSNSIMFLUX_H_INCLUDED__

#include <vector>
#include <utility>
#include <memory>
#include <iostream>
#include <string>

class SKSNSimFluxModel {
  public:
    virtual ~SKSNSimFluxModel() {}
    virtual double /* /cm^2/s */ GetFlux(const double /* MeV */) const = 0; // energy -> flux

    virtual double /* MeV */ GetEnergyLimitMax() const = 0;
    virtual double /* MeV */ GetEnergyLimitMin() const = 0;
    virtual double /* MeV */ GetEnergyLimit(const bool b) const {return ( b? GetEnergyLimitMax(): GetEnergyLimitMin());};
};

class SKSNSimFluxCustom : SKSNSimFluxModel {
  private:
    std::unique_ptr<std::vector<std::pair<double,double>>> ene_flux_v; /* size_t -> <energy, flux> */
    double lower_energy_bin_width;
    double upper_energy_bin_width;
    void sortByEnergy();

    void loadFile(const std::string /*fname*/);
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

  public:
    SKSNSimFluxCustom();
    SKSNSimFluxCustom(const std::string fname);
    ~SKSNSimFluxCustom(){}
    double GetFlux(const double) const;
    double GetEnergyLimitMax() const { return getFluxLimit(false); }
    double GetEnergyLimitMin() const { return getFluxLimit(true); }
    double GetEnergyLimit(const bool b) const { return getFluxLimit(!b); }

};

class SKSNSimFluxSNNakazato : SKSNSimFluxModel {
  public:
    SKSNSimFluxSNNakazato(){}
    ~SKSNSimFluxSNNakazato(){}
    double GetFlux(double e) { return 0.0; }
    double GetEnergyLimitMax() { return 100.0; }
    double GetEnergyLimitMin() { return 10.0; }
};

class SKSNSimFluxDSNBHoriuchi : SKSNSimFluxModel {
  private:
    std::unique_ptr<SKSNSimFluxCustom> customflux;
  public:
    SKSNSimFluxDSNBHoriuchi(){ customflux = std::make_unique<SKSNSimFluxCustom>("./dsnb_flux/horiuchi/8MeV_Nominal.dat"); }
    ~SKSNSimFluxDSNBHoriuchi (){}
    double GetFlux(const double e) const { return customflux->GetFlux(e); }
    double GetEnergyLimitMax() const { return customflux->GetEnergyLimitMax(); }
    double GetEnergyLimitMin() const { return customflux->GetEnergyLimitMin(); }
};

class SKSNSimFluxFlat : SKSNSimFluxModel {
  public:
    SKSNSimFluxFlat (){}
    ~SKSNSimFluxFlat (){}
    double GetFlux(double e) { return 0.0; }
    double GetEnergyLimitMax() { return 100.0; }
    double GetEnergyLimitMin() { return 10.0; }
};

#endif

