/************************************
 * File: SKSNSimUserConfiguration.hh
 ***********************************/

#ifndef __SKSNSIMUSERCONFIGURATION_H_INCLUDED__
#define __SKSNSIMUSERCONFIGURATION_H_INCLUDED__

#include <cstdlib>
#include <string>
#include <vector>
#include <memory>
#include <TRandom3.h>
#include "skrun.h"

class SKSNSimUserConfiguration{
  private:
    double m_energy_min;
    double m_energy_max;
    double m_time_min;
    double m_time_max;
    bool   m_runtime_normalization;
    size_t m_num_events;
    double m_factor_runtime; 
    std::string m_output_directory;
    size_t m_num_per_file;
    std::string m_outputfile_prefix;
    unsigned m_random_seed;
    int m_runtime_runbegin;
    int m_runtime_runend;

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
      m_energy_min = 5.0;
      m_energy_max = 300.0;
      m_time_min = 0.0;
      m_time_max = 0.0;
      m_runtime_normalization = false;
      m_num_events = 1000;
      m_factor_runtime = 1.0;
      m_num_per_file = 1000;
      m_output_directory = "./vectout";
      m_random_seed = 42;
      m_outputfile_prefix = "snmcvector";
      m_runtime_runbegin = SK_VI_BEGIN;
      m_runtime_runend = SK_VI_END;
    }

    SKSNSimUserConfiguration(){
      m_randomgenerator = std::make_shared<TRandom3>();// new TRandom3());
      SetDefaultConfiguation();
    }
    ~SKSNSimUserConfiguration(){}

    void LoadFromArgs(int, char *[]);

    SKSNSimUserConfiguration &SetFluxEnergyMin(double e){ m_energy_min = e; return *this;}
    SKSNSimUserConfiguration &SetFluxEnergyMax(double e){ m_energy_max = e; return *this;}
    SKSNSimUserConfiguration &SetFluxTimeMin(double t){ m_time_min = t; return *this;}
    SKSNSimUserConfiguration &SetFluxTimeMax(double t){ m_time_max = t; return *this;}
    SKSNSimUserConfiguration &SetNumEvents(size_t n){m_num_events = n; return *this;}
    SKSNSimUserConfiguration &SetNormRuntime(bool t){m_runtime_normalization = t; return *this;}
    SKSNSimUserConfiguration &SetRuntimeFactor(double f){ m_factor_runtime = f; return *this;}
    SKSNSimUserConfiguration &SetNumEventsPerFile(size_t n){m_num_per_file = n; return *this;}
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
