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
#include "SKSNSimTools.hh"

class SKSNSimFluxModel {
  public:
    enum FLUXNUTYPE { FLUXNUE = 0, FLUXNUEB, FLUXNUX, NFLUXNUTYPE};
    virtual ~SKSNSimFluxModel() {
      SKSNSimTools::DumpDebugMessage (" dtor SKSNSimFluxModel");
    }
    virtual double /* /cm^2/s */ GetFlux(const double /* MeV */, const double /* sec */, const FLUXNUTYPE) const = 0; // energy -> flux

    virtual double /* MeV */ GetEnergyLimitMax() const = 0;
    virtual double /* MeV */ GetEnergyLimitMin() const = 0;
    virtual double /* sec */ GetTimeLimitMax() const = 0;
    virtual double /* sec */ GetTimeLimitMin() const = 0;
    virtual double /* MeV */ GetEnergyLimit(const bool b) const {return ( b? GetEnergyLimitMax(): GetEnergyLimitMin());};
    virtual double /* sec */ GetTimeLimit(const bool b) const {return ( b? GetTimeLimitMax(): GetTimeLimitMin());};
    virtual const std::set<FLUXNUTYPE> &GetSupportedNuTypes () const = 0;
};

class SKSNSimBinnedFluxModel : public SKSNSimFluxModel {
  public:
    virtual ~SKSNSimBinnedFluxModel() {}
    virtual int GetNBinsEne() const = 0;
    virtual int GetNBinsTime() const = 0;
    virtual double GetBinWidthEne(int b) const = 0;
    virtual double GetBinWidthTime(int b) const = 0;
};

class SKSNSimDSNBFluxCustom : public SKSNSimBinnedFluxModel {
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
    inline double getBinWidth() const {return (ene_flux_v != NULL) ? ene_flux_v->at(1).first - ene_flux_v->at(0).first: -1.0;};
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
    int GetNBinsEne() const { return getNBins(); }
    int GetNBinsTime() const { return 1; }
    double GetBinWidthEne(int b) const { return getBinWidth(); }
    double GetBinWidthTime(int b) const { return 0.0; }
};

class SKSNSimSNFluxCustom : public SKSNSimBinnedFluxModel {
  private:
    std::vector<double> tmesh;
    std::vector<std::vector<double> > enue, eneb, enux, nnue, nneb, nnux, lnue, lneb, lnux;
    const static std::set<FLUXNUTYPE> supportedType;
    int getNBinsEne() const { return enue.front().size(); }
    int getNBinsTime() const { return tmesh.size(); }
    double getBinWidthEne(int b) const { return enue.front().at(1) - enue.front().at(0); }
    double getBinWidthTime(int b) const { return tmesh.at(1) - tmesh.at(0); }
  public:
    SKSNSimSNFluxCustom(){}
    ~SKSNSimSNFluxCustom(){}
    void LoadFluxFile(std::string);
    SKSNSimSNFluxCustom(std::string fname){ LoadFluxFile(fname); }
    double GetFlux(const double e, const double t, const FLUXNUTYPE type) const;
    double GetEnergyLimitMax() const { return enue.front().back(); }
    double GetEnergyLimitMin() const { return enue.front().front(); }
    double GetTimeLimitMax() const { return tmesh.back(); }
    double GetTimeLimitMin() const { return tmesh.front(); }
    const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
    int GetNBinsEne() const { return getNBinsEne(); }
    int GetNBinsTime() const { return getNBinsTime(); }
    double GetBinWidthEne(int b) const { return getBinWidthEne(b); }
    double GetBinWidthTime(int b) const { return getBinWidthTime(b); }
};
const std::set<SKSNSimFluxModel::FLUXNUTYPE> SKSNSimSNFluxCustom::supportedType = {};

class SKSNSimSNFluxNakazatoFormat : public SKSNSimBinnedFluxModel {
  private:
    std::unique_ptr<SKSNSimSNFluxCustom> flux;; // TODO modify to changeable file name (model)

  public:
    SKSNSimSNFluxNakazatoFormat(){
      flux = std::make_unique<SKSNSimSNFluxCustom>( "/home/sklowe/supernova/data/nakazato/intp2002.data" );
    }
    SKSNSimSNFluxNakazatoFormat(std::string mname){
      flux = std::make_unique<SKSNSimSNFluxCustom>( "/home/sklowe/supernova/data/" + mname );
    }
    ~SKSNSimSNFluxNakazatoFormat(){}
    void SetModel(std::string mname) {
      flux.reset(new SKSNSimSNFluxCustom("/home/sklowe/supernova/data/" + mname));
    }
    double GetFlux(const double e, const double t, const FLUXNUTYPE type) const { return flux->GetFlux(e,t,type); }
    double GetEnergyLimitMax() const { return flux->GetEnergyLimitMax(); }
    double GetEnergyLimitMin() const { return flux->GetEnergyLimitMin(); }
    double GetTimeLimitMax() const { return flux->GetTimeLimitMax(); }
    double GetTimeLimitMin() const { return flux->GetTimeLimitMin(); }
    const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return flux->GetSupportedNuTypes(); }
    int GetNBinsEne()        const { return flux->GetNBinsEne(); }
    int GetNBinsTime()       const { return flux->GetNBinsTime(); }
    double GetBinWidthEne(int b)  const { return flux->GetBinWidthEne(b); }
    double GetBinWidthTime(int b) const { return flux->GetBinWidthTime(b); }
    
};

class SKSNSimSNFluxNakazato : public SKSNSimBinnedFluxModel {
  private:
    std::unique_ptr<SKSNSimSNFluxCustom> flux;; // TODO modify to changeable file name (model)
  public:
    SKSNSimSNFluxNakazato(){
      flux = std::make_unique<SKSNSimSNFluxCustom>( "/home/sklowe/supernova/data/nakazato/intp2002.data" );
    }
    ~SKSNSimSNFluxNakazato(){}
    double GetFlux(const double e, const double t, const FLUXNUTYPE type) const { return flux->GetFlux(e,t,type); }
    double GetEnergyLimitMax() const { return flux->GetEnergyLimitMax(); }
    double GetEnergyLimitMin() const { return flux->GetEnergyLimitMin(); }
    double GetTimeLimitMax() const { return flux->GetTimeLimitMax(); }
    double GetTimeLimitMin() const { return flux->GetTimeLimitMin(); }
    const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return flux->GetSupportedNuTypes(); }
    int GetNBinsEne()        const { return flux->GetNBinsEne(); }
    int GetNBinsTime()       const { return flux->GetNBinsTime(); }
    double GetBinWidthEne(int b)  const { return flux->GetBinWidthEne(b); }
    double GetBinWidthTime(int b) const { return flux->GetBinWidthTime(b); }
    
};

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
    double GetEnergyLimitMax() const { return 100.0; }
    double GetEnergyLimitMin() const { return 10.0; }
    const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
};

class SKSNSimFluxMonthlyCustom : SKSNSimFluxModel {
  private: 
    std::vector< std::pair< int /* begin_elapsday */, std::unique_ptr<SKSNSimDSNBFluxCustom> > > custommonthlyflux;
    void sortByTime ();
    const std::set<FLUXNUTYPE> supportedType = {FLUXNUEB};
    SKSNSimDSNBFluxCustom &findFluxByTime(const int /* elapsed day from 1996/01/01 */) const;

  public:
    SKSNSimFluxMonthlyCustom() {}
    ~SKSNSimFluxMonthlyCustom() {}
    void AddMonthlyFlux( const int /* elapse_day from 1996/01/01 */, std::unique_ptr<SKSNSimDSNBFluxCustom> ); /* unique_ptr will be moved to this class */
    double GetFlux(const double /* MeV */, const double /* elapsed day from 1996/01/01 */, const FLUXNUTYPE ) const;
    /* this return flux which fulfilling "t >= (elem[n]->begin_elapsed_day) && t < (elem[n+1]->begin_elapsed_day)" */

    double GetEnergyLimitMax() const { if(custommonthlyflux.empty()) return -1.0;
      return custommonthlyflux.front().second->GetEnergyLimitMax(); }
    double GetEnergyLimitMin() const { if(custommonthlyflux.empty()) return -1.0;
      return custommonthlyflux.front().second->GetEnergyLimitMin(); }
    double GetTimeLimitMax() const { if(custommonthlyflux.empty()) return -1.0;
      return custommonthlyflux.front().first; }
    double GetTimeLimitMin() const { if(custommonthlyflux.empty()) return -1.0;
      return custommonthlyflux.back().first; } // !!Be careful!! THIS IS AN EXCEPTION OF LIMIT VALUE. THIS IS INCLUSIVE.
    const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
};

#endif

