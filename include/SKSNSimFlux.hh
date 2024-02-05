/****************************************
 * File: SKSNSimFlux.hh
 * Description:
 *   Flux model for supernovae simulation
 ****************************************/


#ifndef SKSNSIMFLUX_H_INCLUDED
#define SKSNSIMFLUX_H_INCLUDED

#include <vector>
#include <utility>
#include <memory>
#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include "SKSNSimTools.hh"

constexpr char DATADIRVARIABLENAME[] = "SKSNSIMDATADIR";

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
    virtual double FindMaxFluxTime() const = 0;
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

    void loadFile(const std::string /*fname*/, const std::string /* delimeter */ = "\t");
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
    SKSNSimDSNBFluxCustom(const std::string fname, const std::string /* delimeter */ = "\t");
    ~SKSNSimDSNBFluxCustom(){}
    void DumpFlux(std::ostream &out = std::cout) const;
    double GetFlux(const double, const double t = 0.0, const FLUXNUTYPE nutype = FLUXNUEB) const;
    inline double GetEnergyLimitMax() const { return getFluxLimit(false); }
    inline double GetEnergyLimitMin() const { return getFluxLimit(true); }
    inline double GetEnergyLimit(const bool b) const { return getFluxLimit(!b); }
    inline double GetTimeLimitMax() const { return 0.0;}
    inline double GetTimeLimitMin() const { return 0.0;}
    inline const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
    inline void AddSupoortedNuTypes(std::set<FLUXNUTYPE> s) { for(auto it = s.begin(); it != s.end(); it++) supportedType.emplace(*it);}
    int GetNBinsEne() const { return getNBins(); }
    int GetNBinsTime() const { return 1; }
    double GetBinWidthEne(int b) const { return getBinWidth(); }
    double GetBinWidthTime(int b) const { return 0.0; }
    double CalcIntegratedFlux() const;
    double FindMaxFluxTime() const {return 0.0;}
};

class SKSNSimSNFluxCustom : public SKSNSimBinnedFluxModel {
  private:
    std::vector<double> tmesh;
    std::vector<std::vector<double> > enue, eneb, enux, nnue, nneb, nnux, lnue, lneb, lnux;
    const static std::set<FLUXNUTYPE> supportedType;
    inline int getNBinsEne() const { return enue.front().size(); }
    inline int getNBinsTime() const { return tmesh.size(); }
    inline double getBinWidthEne(int b) const { return enue.front().at(1) - enue.front().at(0); }
    inline double getBinWidthTime(int b) const { return tmesh.at(1) - tmesh.at(0); }
  public:
    SKSNSimSNFluxCustom();
    ~SKSNSimSNFluxCustom(){}
    void LoadFluxFile(std::string);
    SKSNSimSNFluxCustom(std::string fname){ LoadFluxFile(fname); }
    double GetFlux(const double e, const double t, const FLUXNUTYPE type) const;
    inline double GetEnergyLimitMax() const { return enue.front().back(); }
    inline double GetEnergyLimitMin() const { return enue.front().front(); }
    inline double GetTimeLimitMax() const { return tmesh.back(); }
    inline double GetTimeLimitMin() const { return tmesh.front(); }
    inline const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
    inline int GetNBinsEne() const { return getNBinsEne(); }
    inline int GetNBinsTime() const { return getNBinsTime(); }
    inline double GetBinWidthEne(int b) const { return getBinWidthEne(b); }
    inline double GetBinWidthTime(int b) const { return getBinWidthTime(b); } 
    inline double FindMaxFluxTime() const {return 0.0;} // TODO at this momenent, this function does NOT work
};

class SKSNSimSNFluxNakazatoFormat : public SKSNSimBinnedFluxModel {
  private:
    std::unique_ptr<SKSNSimSNFluxCustom> flux;; // TODO modify to changeable file name (model)

  public:
    SKSNSimSNFluxNakazatoFormat(){
      if( const char * env_p = std::getenv(DATADIRVARIABLENAME) )
        flux = std::make_unique<SKSNSimSNFluxCustom>( std::string(env_p) + "/snburst/nakazato/intp2002.data" );
      else {
        std::cout << "The environmental variable \"" << DATADIRVARIABLENAME << "\" is not defined. Please set it..." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    SKSNSimSNFluxNakazatoFormat(std::string mname){
      if( const char * env_p = std::getenv(DATADIRVARIABLENAME) )
        flux = std::make_unique<SKSNSimSNFluxCustom>( std::string(env_p) + "/snburst/" + mname );
      else {
        std::cout << "The environmental variable \"" << DATADIRVARIABLENAME << "\" is not defined. Please set it..." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    ~SKSNSimSNFluxNakazatoFormat(){}
    inline void SetModel(std::string mname) {
      if( const char * env_p = std::getenv(DATADIRVARIABLENAME) )
        flux.reset(new SKSNSimSNFluxCustom( std::string(env_p) + "/snburst/" + mname));
      else {
        std::cout << "The environmental variable \"" << DATADIRVARIABLENAME << "\" is not defined. Please set it..." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    inline double GetFlux(const double e, const double t, const FLUXNUTYPE type) const { return flux->GetFlux(e,t,type); }
    inline double GetEnergyLimitMax() const { return flux->GetEnergyLimitMax(); }
    inline double GetEnergyLimitMin() const { return flux->GetEnergyLimitMin(); }
    inline double GetTimeLimitMax() const { return flux->GetTimeLimitMax(); }
    inline double GetTimeLimitMin() const { return flux->GetTimeLimitMin(); }
    inline const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return flux->GetSupportedNuTypes(); }
    inline int GetNBinsEne()        const { return flux->GetNBinsEne(); }
    inline int GetNBinsTime()       const { return flux->GetNBinsTime(); }
    inline double GetBinWidthEne(int b)  const { return flux->GetBinWidthEne(b); }
    inline double GetBinWidthTime(int b) const { return flux->GetBinWidthTime(b); }
    inline double FindMaxFluxTime() const {return flux->FindMaxFluxTime();}
};

class SKSNSimSNFluxNakazato : public SKSNSimBinnedFluxModel {
  private:
    std::unique_ptr<SKSNSimSNFluxCustom> flux;; // TODO modify to changeable file name (model)
  public:
    SKSNSimSNFluxNakazato(){
      if( const char * env_p = std::getenv(DATADIRVARIABLENAME) )
        flux = std::make_unique<SKSNSimSNFluxCustom>( std::string(env_p) + "/snburst/nakazato/intp2002.data" );
      else {
        std::cout << "The environmental variable \"" << DATADIRVARIABLENAME << "\" is not defined. Please set it..." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    ~SKSNSimSNFluxNakazato(){}
    inline double GetFlux(const double e, const double t, const FLUXNUTYPE type) const { return flux->GetFlux(e,t,type); }
    inline double GetEnergyLimitMax() const { return flux->GetEnergyLimitMax(); }
    inline double GetEnergyLimitMin() const { return flux->GetEnergyLimitMin(); }
    inline double GetTimeLimitMax() const { return flux->GetTimeLimitMax(); }
    inline double GetTimeLimitMin() const { return flux->GetTimeLimitMin(); }
    inline const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return flux->GetSupportedNuTypes(); }
    inline int GetNBinsEne()        const { return flux->GetNBinsEne(); }
    inline int GetNBinsTime()       const { return flux->GetNBinsTime(); }
    inline double GetBinWidthEne(int b)  const { return flux->GetBinWidthEne(b); }
    inline double GetBinWidthTime(int b) const { return flux->GetBinWidthTime(b); }
    inline double FindMaxFluxTime() const {return flux->FindMaxFluxTime();}
};

class SKSNSimFluxDSNBHoriuchi : SKSNSimFluxModel {
  private:
    std::unique_ptr<SKSNSimDSNBFluxCustom> customflux;
  public:
    SKSNSimFluxDSNBHoriuchi(){
      if( const char * env_p = std::getenv(DATADIRVARIABLENAME) )
        customflux = std::make_unique<SKSNSimDSNBFluxCustom>( std::string(env_p) +"/dsnb/horiuchi/8MeV_Nominal.dat");
      else {
        std::cout << "The environmental variable \"" << DATADIRVARIABLENAME << "\" is not defined. Please set it..." << std::endl;
        exit(EXIT_FAILURE);
      }
      customflux->AddSupoortedNuTypes(std::set<FLUXNUTYPE>({FLUXNUEB}));
    }
    ~SKSNSimFluxDSNBHoriuchi (){}
    inline double GetFlux(const double e, const double t, const FLUXNUTYPE type) const { return customflux->GetFlux(e,t, type); }
    inline double GetEnergyLimitMax() const { return customflux->GetEnergyLimitMax(); }
    inline double GetEnergyLimitMin() const { return customflux->GetEnergyLimitMin(); }
    inline double GetTimeLimitMax() const { return customflux->GetEnergyLimitMax(); }
    inline double GetTimeLimitMin() const { return customflux->GetEnergyLimitMin(); }
    inline void DumpFlux(std::ostream &out = std::cout) const { customflux->DumpFlux(out);};
    inline const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return customflux->GetSupportedNuTypes(); }
    inline double FindMaxFluxTime() const {return customflux->FindMaxFluxTime();}
};

class SKSNSimFluxFlat : SKSNSimFluxModel {
  private:
    std::set<FLUXNUTYPE> supportedType;
  public:
    SKSNSimFluxFlat (){}
    ~SKSNSimFluxFlat (){}
    inline double GetFlux(const double e, const double t, const FLUXNUTYPE type) const { return 0.0; }
    inline double GetEnergyLimitMax() const { return 100.0; }
    inline double GetEnergyLimitMin() const { return 10.0; }
    inline const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
    inline double FindMaxFluxTime() const {return 0.0;}
};

class SKSNSimDSNBFluxMonthlyCustom : SKSNSimFluxModel {
  private: 
    std::vector< std::pair< int /* begin_elapsday */, std::unique_ptr<SKSNSimDSNBFluxCustom> > > custommonthlyflux;
    void sortByTime ();
    const static std::set<FLUXNUTYPE> supportedType;
    SKSNSimDSNBFluxCustom &findFluxByTime(const int /* elapsed day from 1996/01/01 */) const;

  public:
    SKSNSimDSNBFluxMonthlyCustom();
    ~SKSNSimDSNBFluxMonthlyCustom() {}
    void AddMonthlyFlux( const int /* elapse_day from 1996/01/01 */, std::unique_ptr<SKSNSimDSNBFluxCustom> ); /* unique_ptr will be moved to this class */
    double GetFlux(const double /* MeV */, const double /* elapsed day from 1996/01/01 */, const FLUXNUTYPE ) const;
    /* this return flux which fulfilling "t >= (elem[n]->begin_elapsed_day) && t < (elem[n+1]->begin_elapsed_day)" */

    double FindMaxFluxTime() const ;
    inline SKSNSimDSNBFluxCustom &FindFluxByTime(const int d /* elapsed day from 1996/01/01 */) const { return findFluxByTime(d);}

    inline double GetEnergyLimitMax() const { if(custommonthlyflux.empty()) return -1.0;
      return custommonthlyflux.front().second->GetEnergyLimitMax(); }
    inline double GetEnergyLimitMin() const { if(custommonthlyflux.empty()) return -1.0;
      return custommonthlyflux.front().second->GetEnergyLimitMin(); }
    inline double GetTimeLimitMax() const { if(custommonthlyflux.empty()) return -1.0;
      return custommonthlyflux.front().first; }
    inline double GetTimeLimitMin() const { if(custommonthlyflux.empty()) return -1.0;
      return custommonthlyflux.back().first; } // !!Be careful!! THIS IS AN EXCEPTION OF LIMIT VALUE. THIS IS INCLUSIVE.
    inline const std::set<FLUXNUTYPE> &GetSupportedNuTypes() const { return supportedType; }
};

#endif

