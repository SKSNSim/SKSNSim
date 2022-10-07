/************************************
 * File: SKSNSimUserConfiguration.hh
 * Author: Shota Izumiyama (izumiyama@hep.phys.titech.ac.jp)
 * Date: Fri Oct  7 10:08:11 JST 2022
 * Desctiption:
 *   Configuration
 ***********************************/

#ifndef __SKSNSIMUSERCONFIGURATION_H_INCLUDED__
#define __SKSNSIMUSERCONFIGURATION_H_INCLUDED__

#include <cstdlib>
#include <string>
#include <vector>
#include <memory>
#include <TRandom3.h>
#include "skrun.h"
#include "SKSNSimEnum.hh"

class SKSNSimUserConfiguration{
  private:

    /* Event range related */
    double m_energy_min;
    double m_energy_max;
    double m_time_min;
    double m_time_max;

    /* Output file related */
    std::string m_output_directory;
    std::string m_outputfile_prefix;
    size_t m_num_per_file;
    bool m_eventvector_generation;
    SKSNSIMENUM::TANKVOLUME m_eventgen_volume;

    /* DSNB runtime normalization related */
    size_t m_num_events;
    bool   m_runtime_normalization;
    double m_factor_runtime; 
    int m_runtime_runbegin;
    int m_runtime_runend;

    /* Physics related */
    double m_sndistance_kpc;
    SKSNSIMENUM::NEUTRINOOSCILLATION m_nuosc_type;
    std::string m_snburst_fluxmodel;

    /* Random Generator related */
    unsigned m_random_seed;
    std::shared_ptr<TRandom> m_randomgenerator;

    bool CheckFluxEnergyMin() const;
    bool CheckFluxEnergyMax() const;
    bool CheckFluxTimeMin() const;
    bool CheckFluxTimeMax() const;
    bool CheckNumEvents() const;
    bool CheckNormRuntime() const;
    bool CheckRuntimeFactor() const;

  public:
    void SetDefaultConfiguation() {
      m_energy_min = GetDefaultFluxEnergyMin();
      m_energy_max = GetDefaultFluxEnergyMax();
      m_time_min = GetDefaultFluxTimeMin();
      m_time_max = GetDefaultFluxTimeMax();

      m_output_directory = GetDefaultOutputDirectory();
      m_outputfile_prefix = GetDefaultOutputPrefix();
      m_num_per_file = GetDefaultNumEventsPerFile();
      m_eventvector_generation = GetDefaultVectorGeneration();
      m_eventgen_volume = GetDefaultEventVolume();

      m_num_events = GetDefaultNumEvents();
      m_runtime_normalization = GetDefaultRuttimeNormalization();
      m_factor_runtime = GetDefaultRuttimeNormFactor();
      m_runtime_runbegin = GetDefaultRuntimeBegin();
      m_runtime_runend = GetDefaultRuntimeEnd();

      m_sndistance_kpc = GetDefaultSNDistanceKPC();
      m_nuosc_type = GetDefaultNeutrinoOscType();
      m_snburst_fluxmodel = GetDefaultSNBurstFluxModel();

      m_random_seed = GetDefaultRandomSeed();
    }

    SKSNSimUserConfiguration(){
      m_randomgenerator = std::make_shared<TRandom3>();// new TRandom3());
      SetDefaultConfiguation();
    }
    ~SKSNSimUserConfiguration(){}

    void LoadFromArgs(int, char *[]); // Parser from arguments of the main function

    /* Default configuration: please define default values through them */
    const static SKSNSIMENUM::SIMTYPE GetDefaultSimType() { return SKSNSIMENUM::SIMTYPE::kSNBURST; }
    const static double GetDefaultFluxEnergyMin () { return   5.0 /* MeV */;}
    const static double GetDefaultFluxEnergyMax () { return 300.0 /* MeV */;}
    const static double GetDefaultFluxTimeMin () { return   0. /* sec */;}
    const static double GetDefaultFluxTimeMax () { return 200. /* sec */;}
    const static SKSNSIMENUM::TANKVOLUME GetDefaultEventVolume () { return SKSNSIMENUM::TANKVOLUME::kIDFV; }
    const static SKSNSIMENUM::NEUTRINOOSCILLATION GetDefaultNeutrinoOscType () { return SKSNSIMENUM::NEUTRINOOSCILLATION::kNONE; }
    const static double GetDefaultSNDistanceKPC () { return 10. /* kpc */;}
    const static bool GetDefaultVectorGeneration () { return true; } 
    const static std::string GetDefaultOutputDirectory () { return "./vectout"; }
    const static std::string GetDefaultOutputPrefix () { return "snmcvect"; }
    const static size_t GetDefaultNumEvents () { return 1000;}
    const static size_t GetDefaultNumEventsPerFile () { return 1000; }
    const static unsigned GetDefaultRandomSeed () { return 42;}
    const static int GetDefaultRuntimeBegin () { return (int) SKSNSIMENUM::SKPERIODRUN::SKVIBEGIN;}
    const static int GetDefaultRuntimeEnd () { return (int) SKSNSIMENUM::SKPERIODRUN::SKVIEND;}
    const static bool GetDefaultRuttimeNormalization() { return false; }
    const static double GetDefaultRuttimeNormFactor() { return 1.0; }
    const static std::string GetDefaultSNBurstFluxModel () { return "nakazato/intp2002.data";}


    /* =========================
     * Setting methods
     * If they returns myselfs, we can write the codes like this:
     *   config.SetFluxEnergyMin(0.0)
     *         .SetFluxEnergyMax(300.0);
     * ==========================*/
    SKSNSimUserConfiguration &SetFluxEnergyMin(double e){ m_energy_min = e; return *this;}
    SKSNSimUserConfiguration &SetFluxEnergyMax(double e){ m_energy_max = e; return *this;}
    SKSNSimUserConfiguration &SetFluxTimeMin(double t){ m_time_min = t; return *this;}
    SKSNSimUserConfiguration &SetFluxTimeMax(double t){ m_time_max = t; return *this;}
    SKSNSimUserConfiguration &SetNumEvents(size_t n){m_num_events = n; return *this;}
    SKSNSimUserConfiguration &SetNumEventsPerFile(size_t n){m_num_per_file = n; return *this;}
    SKSNSimUserConfiguration &SetNormRuntime(bool t){m_runtime_normalization = t; return *this;}
    SKSNSimUserConfiguration &SetRuntimeFactor(double f){ m_factor_runtime = f; return *this;}
    SKSNSimUserConfiguration &SetOutputDirectory(std::string dir) { m_output_directory = dir; return *this;}
    SKSNSimUserConfiguration &SetOutputPrefix(std::string pref) { m_outputfile_prefix = pref; return *this;}
    SKSNSimUserConfiguration &SetRandomSeed(unsigned s) { m_random_seed = s; m_randomgenerator->SetSeed(m_random_seed); return *this;}
    SKSNSimUserConfiguration &SetRuntimeRunBegin(int r) { m_runtime_runbegin = r; return *this;}
    SKSNSimUserConfiguration &SetRuntimeRunEnd(int r) { m_runtime_runend = r; return *this;}


    double GetFluxEnergyMin() const { return m_energy_min;}
    double GetFluxEnergyMax() const { return m_energy_max;}
    double GetFluxTimeMin() const { return m_time_min;}
    double GetFluxTimeMax() const { return m_time_max;}
    bool   GetNormRuntime() const { return m_runtime_normalization;}
    size_t GetNumEvents() const { return m_num_events;}
    double GetRuntimeFactor() const { return m_factor_runtime;}
    size_t GetNumEventsPerFile() const {return m_num_per_file;}
    unsigned GetRandomSeed() const {return m_random_seed;}
    std::string GetOutputDirectory() const { return m_output_directory;}
    std::string GetOutputPrefix() const { return m_outputfile_prefix;}
    std::shared_ptr<TRandom> GetRandomGenerator() { return m_randomgenerator;}
    int GetRuntimeRunBegin() const {return m_runtime_runbegin; }
    int GetRuntimeRunEnd() const {return m_runtime_runend; }

    bool CheckHealth () const;
    void Dump() const;
};


#endif
