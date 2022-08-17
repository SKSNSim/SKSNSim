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
#include <set>

class SKSNSimFluxModel {
  public:
    enum FLUXNUTYPE { FLUXNUE = 0, FLUXNUEB, FLUXNUX, NFLUXNUTYPE};
    virtual ~SKSNSimFluxModel() {}
    virtual double /* /cm^2/s */ GetFlux(const double /* MeV */, const double /* sec */, const FLUXNUTYPE) const = 0; // energy -> flux

    virtual double /* MeV */ GetEnergyLimitMax() const = 0;
    virtual double /* MeV */ GetEnergyLimitMin() const = 0;
    virtual double /* sec */ GetTimeLimitMax() const = 0;
    virtual double /* sec */ GetTimeLimitMin() const = 0;
    virtual double /* MeV */ GetEnergyLimit(const bool b) const {return ( b? GetEnergyLimitMax(): GetEnergyLimitMin());};
    virtual double /* sec */ GetTimeLimit(const bool b) const {return ( b? GetTimeLimitMax(): GetTimeLimitMin());};
    virtual const std::set<FLUXNUTYPE> &GetSupportedNuTypes () const = 0;
};

class SKSNSimDSNBFluxCustom : SKSNSimFluxModel {
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
    inline double getBinnedFlux(const int bin) const {return ene_flux_v->at(bin).second;};
    inline double getBinnedEnergy(const int bin) const {return ene_flux_v->at(bin).first;};
    inline size_t getNBins() const {return (ene_flux_v != NULL) ? ene_flux_v->size(): 0;};
    std::set<FLUXNUTYPE> supportedType;

  public:
    SKSNSimDSNBFluxCustom();
    SKSNSimDSNBFluxCustom(const std::string fname);
    ~SKSNSimDSNBFluxCustom(){}
    void DumpFlux(std::ostream &out = std::cout) const;
    double GetFlux(const double, const double t = 0.0, const FLUXNUTYPE nutype = FLUXNUEB) const;
    double GetEnergyLimitMax() const { return getFluxLimit(false); }
    double GetEnergyLimitMin() const { return getFluxLimit(true); }
    double GetEnergyLimit(const bool b) const { return getFluxLimit(!b); }
    double GetTimeLimitMax() const { return 0.0;}
    double GetTimeLimitMin() const { return 0.0;}
    const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
    void AddSupoortedNuTypes(std::set<FLUXNUTYPE> s) { for(auto it = s.begin(); it != s.end(); it++) supportedType.emplace(*it);}

};

class SKSNSimSNFluxNakazato : SKSNSimFluxModel {
  private:
    const static std::set<FLUXNUTYPE> supportedType;
  public:
    SKSNSimSNFluxNakazato(){}
    ~SKSNSimSNFluxNakazato(){}
    double GetFlux(const double e, const double t, const FLUXNUTYPE type) const { return 0.0; }
    double GetEnergyLimitMax() { return 100.0; }
    double GetEnergyLimitMin() { return 10.0; }
    double GetTimeLimitMax() { return 100.0; }
    double GetTimeLimitMin() { return 0.0; }
    const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
};
const std::set<SKSNSimFluxModel::FLUXNUTYPE> SKSNSimSNFluxNakazato::supportedType = { SKSNSimFluxModel::FLUXNUEB, SKSNSimFluxModel::FLUXNUE, SKSNSimFluxModel::FLUXNUX };

class SKSNSimFluxDSNBHoriuchi : SKSNSimFluxModel {
  private:
    std::unique_ptr<SKSNSimDSNBFluxCustom> customflux;
  public:
    SKSNSimFluxDSNBHoriuchi(){
      customflux = std::make_unique<SKSNSimDSNBFluxCustom>("./dsnb_flux/horiuchi/8MeV_Nominal.dat");
      customflux->AddSupoortedNuTypes(std::set<FLUXNUTYPE>({FLUXNUEB}));
    }
    ~SKSNSimFluxDSNBHoriuchi (){}
    double GetFlux(const double e, const double t, const FLUXNUTYPE type) const { return customflux->GetFlux(e,t, type); }
    double GetEnergyLimitMax() const { return customflux->GetEnergyLimitMax(); }
    double GetEnergyLimitMin() const { return customflux->GetEnergyLimitMin(); }
    double GetTimeLimitMax() const { return customflux->GetEnergyLimitMax(); }
    double GetTimeLimitMin() const { return customflux->GetEnergyLimitMin(); }
    void DumpFlux(std::ostream &out = std::cout) const { customflux->DumpFlux(out);};
    const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return customflux->GetSupportedNuTypes(); }
};

class SKSNSimFluxFlat : SKSNSimFluxModel {
  private:
    std::set<FLUXNUTYPE> supportedType;
  public:
    SKSNSimFluxFlat (){}
    ~SKSNSimFluxFlat (){}
    double GetFlux(const double e, const double t, const FLUXNUTYPE type) const { return 0.0; }
    double GetEnergyLimitMax() { return 100.0; }
    double GetEnergyLimitMin() { return 10.0; }
    const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
};

#endif

