/************************************
 * File: SKSNSimUserConfiguration.hh
 * Author: Shota Izumiyama (izumiyama@hep.phys.titech.ac.jp)
 * Date: Fri Oct  7 10:08:11 JST 2022
 * Desctiption:
 *   Configuration
 ***********************************/

#ifndef SKSNSIMUSERCONFIGURATION_H_INCLUDED
#define SKSNSIMUSERCONFIGURATION_H_INCLUDED

#include <cstdlib>
#include <string>
#include <vector>
#include <memory>
#include <TRandom3.h>
#include "skrun.h"
#include "SKSNSimEnum.hh"
#include "SKSNSimVectorGenerator.hh"


class SKSNSimUserConfiguration{
  public:
    enum struct MODERUNTIME   { kEVNUM = 0, kRUNTIMERUNNUM, kRUNTIMEPERIOD, kNMODERUNTIME };
    enum struct MODEGENERATOR { kSNBURST = 0, kDSNB, kNMODEGENERATOR };
    enum struct MODEOFILE { kSKROOT = 0, kNUANCE, kNMODEOFILE };

  private:
    /* mode */
    const MODEGENERATOR m_mode_generator;

    /* Event range related */
    double m_energy_min;
    double m_energy_max;
    size_t m_energy_nbins;
    double m_time_min;
    double m_time_max;
    size_t m_time_nbins;

    /* Output file related */
    std::string m_output_directory;
    std::string m_outputfile_prefix;
    std::string m_outputfile_template;
    size_t m_num_per_file;
    bool m_eventvector_generation;
    SKSNSIMENUM::TANKVOLUME m_eventgen_volume;
    MODEOFILE m_mode_ofile;

    /* DSNB runtime normalization related */
    size_t m_num_events;
    bool   m_runtime_normalization;
    double m_factor_runtime; 
    int m_runtime_runbegin;
    int m_runtime_runend;
    int m_runtime_period;

    /* RUN number (exclusive to DSNB runtime normalization) */
    int m_runnum;
    int m_subrunnum;

    /* SN related */
    double m_sndistance_kpc;

    /* Physics related */
    SKSNSIMENUM::NEUTRINOOSCILLATION m_nuosc_type;
    std::string m_snburst_fluxmodel;
    std::string m_dsnb_fluxmodel;
    bool m_dsnb_flatflux;

    /* Random Generator related */
    unsigned m_random_seed;
    std::shared_ptr<TRandom> m_randomgenerator;

    /* Health checker. They should return true (if configuration is good (false if bad)) */
    bool CheckFluxEnergyMin() const;
    bool CheckFluxEnergyMax() const;
    bool CheckFluxTimeMin() const;
    bool CheckFluxTimeMax() const;
    bool CheckNumEvents() const;
    bool CheckNormRuntime() const;
    bool CheckRuntimeFactor() const;
    bool CheckOFileMode() const { return m_mode_ofile != MODEOFILE::kNMODEOFILE; }

    static std::string convOFileModeString(MODEOFILE m);

  public:
    void SetDefaultConfiguation() {
      m_energy_min = GetDefaultFluxEnergyMin(m_mode_generator);
      m_energy_max = GetDefaultFluxEnergyMax(m_mode_generator);
      m_energy_nbins = GetDefaultEnergyNBins();
      m_time_min = GetDefaultFluxTimeMin();
      m_time_max = GetDefaultFluxTimeMax();
      m_time_nbins = GetDefaultTimeNBins();

      m_output_directory = GetDefaultOutputDirectory();
      m_outputfile_prefix = GetDefaultOutputPrefix();
      m_outputfile_template = GetDefaultOutputNameTemplate();
      m_num_per_file = GetDefaultNumEventsPerFile();
      m_eventvector_generation = GetDefaultVectorGeneration();
      m_eventgen_volume = GetDefaultEventVolume();
      m_mode_ofile = GetDefaultOFileMode();

      m_num_events = GetDefaultNumEvents();
      m_runtime_normalization = GetDefaultRuntimeNormalization();
      m_factor_runtime = GetDefaultRuntimeNormFactor();
      m_runtime_runbegin = GetDefaultRuntimeBegin();
      m_runtime_runend = GetDefaultRuntimeEnd();
      m_runtime_period = GetDefaultRuntimePeriod();

      m_runnum = GetDefaultRunnum();
      m_subrunnum = GetDefaultSubRunnum();

      m_sndistance_kpc = GetDefaultSNDistanceKPC();
      m_snburst_fluxmodel = GetDefaultSNBurstFluxModel();
      m_dsnb_fluxmodel = GetDefaultDSNBFluxModel();
      m_dsnb_flatflux = GetDefaultDSNBFlatFlux();

      m_nuosc_type = GetDefaultNeutrinoOscType();

      m_random_seed = GetDefaultRandomSeed();
    }

    SKSNSimUserConfiguration( MODEGENERATOR mode = MODEGENERATOR::kDSNB ): 
      m_mode_generator( mode ) {
        m_randomgenerator = std::make_shared<TRandom3>( GetDefaultRandomSeed() );// new TRandom3());
        SetDefaultConfiguation();
      }
    ~SKSNSimUserConfiguration(){}

    void LoadFromArgsDSNB(int, char *[]); // Parser from arguments of the main function. For DSNB  generator program
    void LoadFromArgsSN(int, char *[]); // Parser from arguments of the main function this is for SN burst program

    /* Default configuration: please define default values through them */
    const static double GetDefaultFluxEnergyMin (MODEGENERATOR m) { if(m==MODEGENERATOR::kSNBURST) return   0.0 /* MeV */;
      else if(m==MODEGENERATOR::kDSNB) return 5.0 /* MeV */;
      return 0.0 /* MeV */;
    }
    const static double GetDefaultFluxEnergyMax (MODEGENERATOR m) { if(m==MODEGENERATOR::kSNBURST) return 300.0 /* MeV */;
      else if(m==MODEGENERATOR::kDSNB) return 80.0 /* MeV */;
      return 300.0 /* MeV */;
    }
    const static size_t GetDefaultEnergyNBins () { return 3000;}
    const static double GetDefaultFluxTimeMin () { return   0. /* sec */;}
    const static double GetDefaultFluxTimeMax () { return 20. /* sec */;}
    const static size_t GetDefaultTimeNBins () { return 20000;}
    const static SKSNSIMENUM::TANKVOLUME GetDefaultEventVolume () { return SKSNSIMENUM::TANKVOLUME::kIDFULL; }
    const static SKSNSIMENUM::NEUTRINOOSCILLATION GetDefaultNeutrinoOscType () { return SKSNSIMENUM::NEUTRINOOSCILLATION::kNONE; }
    const static double GetDefaultSNDistanceKPC () { return 10. /* kpc */;}
    const static std::string GetDefaultSNModelName () { return "nakazato/intp2002.data" ;}
    const static bool GetDefaultVectorGeneration () { return true; } 
    const static std::string GetDefaultOutputDirectory () { return "./vectout"; }
    const static std::string GetDefaultOutputPrefix () { return "snmcvect"; }
    const static std::string GetDefaultOutputNameTemplate () { return ""; }
    const static size_t GetDefaultNumEvents () { return 1000;}
    const static size_t GetDefaultNumEventsPerFile () { return 1000; }
    const static unsigned GetDefaultRandomSeed () { return 42;}
    const static MODEOFILE GetDefaultOFileMode () { return MODEOFILE::kSKROOT; }
    const static int GetDefaultRuntimeBegin () { return (int) SKSNSIMENUM::SKPERIODRUN::SKVIBEGIN;}
    const static int GetDefaultRuntimeEnd () { return (int) SKSNSIMENUM::SKPERIODRUN::SKVIEND;}
    const static int GetDefaultRuntimePeriod() { return -1; /* using RuntimeBegin/End */}
    const static bool GetDefaultRuntimeNormalization() { return false; }
    const static double GetDefaultRuntimeNormFactor() { return 24.0; }
    const static std::string GetDefaultSNBurstFluxModel () {
      std::string dir;
      if( const char * env_p = std::getenv(DATADIRVARIABLENAME) )
        dir = std::string(env_p);
      else {
        std::cout << "The environmental variable default flux model\"" << DATADIRVARIABLENAME << "\" is not defined. Please set it..." << std::endl;
        exit(EXIT_FAILURE);
      }
      //return dir + "/nakazato/intp2002.data";
      return "/nakazato/intp2002.data";
    }
    const static std::string GetDefaultDSNBFluxModel () {
      std::string dir;
      if( const char * env_p = std::getenv(DATADIRVARIABLENAME) )
        dir = std::string(env_p);
      else {
        std::cout << "The environmental variable \"" << DATADIRVARIABLENAME << "\" is not defined. Please set it..." << std::endl;
        exit(EXIT_FAILURE);
      }
      return dir + "/horiuchi/8MeV_Nominal.dat";
    }
    const static bool GetDefaultDSNBFlatFlux () { return false;}
    const static int GetDefaultRunnum () { return (int) SKSNSIMENUM::SKPERIODRUN::SKMC; }
    const static int GetDefaultSubRunnum () { return 0; }


    /* =========================
     * Setting methods
     * If they returns myselfs, we can write the codes like this:
     *   config.SetFluxEnergyMin(0.0)
     *         .SetFluxEnergyMax(300.0);
     * ==========================*/
    SKSNSimUserConfiguration &SetFluxEnergyMin(double e){ m_energy_min = e; return *this;}
    SKSNSimUserConfiguration &SetFluxEnergyMax(double e){ m_energy_max = e; return *this;}
    SKSNSimUserConfiguration &SetEnergyNBins(size_t n){ m_energy_nbins = n; return *this;}
    SKSNSimUserConfiguration &SetFluxTimeMin(double t){ m_time_min = t; return *this;}
    SKSNSimUserConfiguration &SetFluxTimeMax(double t){ m_time_max = t; return *this;}
    SKSNSimUserConfiguration &SetTimeNBins(size_t n){ m_time_nbins = n; return *this;}
    SKSNSimUserConfiguration &SetNumEvents(size_t n){m_num_events = n; return *this;}
    SKSNSimUserConfiguration &SetNumEventsPerFile(size_t n){m_num_per_file = n; return *this;}
    SKSNSimUserConfiguration &SetNormRuntime(bool t){m_runtime_normalization = t; return *this;}
    SKSNSimUserConfiguration &SetRuntimeFactor(double f){ m_factor_runtime = f; return *this;}
    SKSNSimUserConfiguration &SetOutputDirectory(std::string dir) { m_output_directory = dir; return *this;}
    SKSNSimUserConfiguration &SetOutputPrefix(std::string pref) { m_outputfile_prefix = pref; return *this;}
    SKSNSimUserConfiguration &SetOutputNameTemplate(std::string tmp) { m_outputfile_template = tmp; return *this;}
    SKSNSimUserConfiguration &SetRandomSeed(unsigned s) { m_random_seed = s; m_randomgenerator->SetSeed(m_random_seed); return *this;}
    SKSNSimUserConfiguration &SetRuntimeRunBegin(int r) { m_runtime_runbegin = r; return *this;}
    SKSNSimUserConfiguration &SetRuntimeRunEnd(int r) { m_runtime_runend = r; return *this;}
    SKSNSimUserConfiguration &SetRuntimePeriod(int p) { m_runtime_period = p; return *this;}
    SKSNSimUserConfiguration &SetSNDistanceKpc(double d) { m_sndistance_kpc = d; return *this;}
    SKSNSimUserConfiguration &SetSNBurstFluxModel(std::string f) { m_snburst_fluxmodel = f; return *this;}
    SKSNSimUserConfiguration &SetDSNBFluxModel(std::string f) { m_dsnb_fluxmodel = f; return *this;}
    SKSNSimUserConfiguration &SetDSNBFlatFlux(bool f) { m_dsnb_flatflux = f; return *this;}
    SKSNSimUserConfiguration &SetVectorGeneration(bool f) { m_eventvector_generation = f; return *this;}
    SKSNSimUserConfiguration &SetNeutrinoOscType( SKSNSIMENUM::NEUTRINOOSCILLATION t) { m_nuosc_type = t; return *this; }
    SKSNSimUserConfiguration &SetNeutrinoOscType( int t) { m_nuosc_type = (SKSNSIMENUM::NEUTRINOOSCILLATION)t; return *this; }
    SKSNSimUserConfiguration &SetRunnum(int r) { m_runnum = r; return *this; }
    SKSNSimUserConfiguration &SetSubRunnum(int r) { m_subrunnum = r; return *this; }
    SKSNSimUserConfiguration &SetOFileMode ( MODEOFILE m ) { m_mode_ofile = m; return *this; }
    SKSNSimUserConfiguration &SetOFileMode ( std::string s, bool exit_if_wrong = true );

    /* Event range related */
    double GetFluxEnergyMin() const { return m_energy_min;}
    double GetFluxEnergyMax() const { return m_energy_max;}
    size_t GetEnergyNBins() const { return m_energy_nbins;}
    double GetFluxTimeMin() const { return m_time_min;}
    double GetFluxTimeMax() const { return m_time_max;}
    size_t GetTimeNBins() const { return m_time_nbins;}

    /* Output file related */
    std::string GetOutputDirectory() const { return m_output_directory;}
    std::string GetOutputPrefix() const { return m_outputfile_prefix;}
    std::string GetOutputNameTemplate() const { return m_outputfile_template;}
    size_t GetNumEventsPerFile() const {return m_num_per_file;}
    bool   GetEventVectorGeneration() const { return m_eventvector_generation; }
    SKSNSIMENUM::TANKVOLUME GetEventgenVolume() const { return m_eventgen_volume; }
    MODEOFILE GetOFileMode() const { return m_mode_ofile; }
    std::string GetOFileModeString() const;

    /* DSNB runtime normalization related */
    bool   GetNormRuntime() const { return m_runtime_normalization;}
    size_t GetNumEvents() const { return m_num_events;}
    double GetRuntimeFactor() const { return m_factor_runtime;}
    int GetRuntimeRunBegin() const {return m_runtime_runbegin; }
    int GetRuntimeRunEnd() const {return m_runtime_runend; }
    int GetRuntimePeriod() const { return m_runtime_period; }

    /* RUN number (exclusive to DSNB runtime normalization) */
    int GetRunnum() const { return m_runnum;}
    int GetSubRunnum() const { return m_subrunnum;}

    /* SN related */
    double GetSNDistanceKpc() const { return m_sndistance_kpc; }
    std::string GetSNBurstFluxModel() const { return m_snburst_fluxmodel; }
    std::string GetDSNBFluxModel() const { return m_dsnb_fluxmodel; }
    bool GetDSNBFlatFlux() const { return m_dsnb_flatflux; }

    /* Physics related */
    SKSNSIMENUM::NEUTRINOOSCILLATION GetNuOscType() const { return m_nuosc_type; }

    unsigned GetRandomSeed() const {return m_random_seed;}
    std::shared_ptr<TRandom> GetRandomGenerator() { return m_randomgenerator;}

    bool CheckHealth () const;
    void Dump() const;


    MODERUNTIME CheckMODERuntime() const;

    void Apply( SKSNSimVectorSNGenerator &gen ) const ;
    void Apply( SKSNSimVectorGenerator   &gen ) const ;

  public:
    static void ShowHelpDSNB(const char *argv0);
    static void ShowHelpSN(const char *argv0);
};


#endif
