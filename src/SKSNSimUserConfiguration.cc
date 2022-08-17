/************************************
 * File: SKSNSimUserConfiguration.cc
 ***********************************/

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include "SKSNSimUserConfiguration.hh"

bool SKSNSimUserConfiguration::CheckFluxEnergyMin() const {
  bool badhealth = false;
  if(m_energy_min > m_energy_max) {
    std::cerr << "FluxEnergyMin: bad (" << m_energy_min << "): smaller than max (" << m_energy_max << ") " << std::endl;
    badhealth |= true;
  }
  if(m_energy_min < 0.){
    std::cerr << "FluxEnergyMin: bad (" << m_energy_min << "): out of range " << std::endl;
    badhealth |= true;
  }
  return !badhealth;
}
bool SKSNSimUserConfiguration::CheckFluxEnergyMax() const {
  bool badhealth = false;
  if(m_energy_min > m_energy_max) {
    std::cerr << "FluxEnergyMax: bad (" << m_energy_max << "): larger than min (" << m_energy_min << ") " << std::endl;
    badhealth |= true;
  }
  if(m_energy_max < 0.){
    std::cerr << "FluxEnergyMin: bad (" << m_energy_max << "): out of range " << std::endl;
    badhealth |= true;
  }
  return !badhealth;
}
bool SKSNSimUserConfiguration::CheckFluxTimeMin() const {
  bool badhealth = false;
  if(m_time_min > m_time_max) {
    std::cerr << "FluxTimeMin: bad (" << m_time_min << "): smaller than max (" << m_time_max << ") " << std::endl;
    badhealth |= true;
  }
  if(m_time_min < 0.){
    std::cerr << "FluxTimeMin: bad (" << m_time_min << "): out of range " << std::endl;
    badhealth |= true;
  }
  return !badhealth;
}
bool SKSNSimUserConfiguration::CheckFluxTimeMax() const {
  bool badhealth = false;
  if(m_time_min > m_time_max) {
    std::cerr << "FluxTimeMax: bad (" << m_time_max << "): larger than min (" << m_time_min << ") " << std::endl;
    badhealth |= true;
  }
  if(m_time_max < 0.){
    std::cerr << "FluxTimeMin: bad (" << m_time_max << "): out of range " << std::endl;
    badhealth |= true;
  }
  return !badhealth;
}
bool SKSNSimUserConfiguration::CheckNumEvents() const {
  bool badhealth = false;
  if(m_num_events <= 0){
    std::cerr << "NumEvents: bad (" << m_num_events << "): == 0" << std::endl;
    badhealth |= true;
  }
  return !badhealth;
}
bool SKSNSimUserConfiguration::CheckNormRuntime() const {
  bool badhealth = false;
  if( m_runtime_normalization && m_factor_runtime <= 0.0 ){
    std::cerr << "NormRuntime: bad (" << m_runtime_normalization << "): factor is not set correctly (factor=" << m_factor_runtime << ")" << std::endl;
    badhealth |= true;
  }
  if( !m_runtime_normalization && m_num_events <= 0){
    std::cerr << "NormRuntime: bad (" << m_runtime_normalization << "): NumEvents is not set correctly (NumEvents=" << m_num_events << ")" << std::endl;
    badhealth |= true;
  }
  return !badhealth;
}
bool SKSNSimUserConfiguration::CheckRuntimeFactor() const {
  bool badhealth = false;
  if( m_factor_runtime <= 0.0 ){
    std::cerr << "RuntimeFactor: bad (" << m_factor_runtime << "): out of range" << std::endl;
    badhealth |= true;
  }
  return !badhealth;
}

void ShowHelp(const char *argv0){
  std::cout << argv0
    << " [-h,--help]"
    << " [-o,--outdir outputdirectory]"
    << " [-n,--nevents numberofevents]"
    << " [--neventsperfile numberofevents]"
    << " [-s,--seed randomseed]"
    << " [--energy_min energy_MeV]"
    << " [--energy_max energy_MeV]"
    << " [--time_min time]"
    << " [--time_max time]"
    << " [--runtimefactor factor]"
    << " [--runtime]"
    << " [--outprefix prefix]"
    << " outputdirectory"
    << std::endl;
}

void SKSNSimUserConfiguration::LoadFromArgs(int argc, char *argv[]){
  int c;
  while(1){
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;

    static struct option long_options[] = {
      {"energy_min",    required_argument, 0,   0},
      {"energy_max",    required_argument, 0,   0},
      {"time_min",      required_argument, 0,   0},
      {"time_max",      required_argument, 0,   0},
      {"runtimefactor", required_argument, 0,   0},
      {"nevents",       required_argument, 0, 'n'},
      {"outdir",        required_argument, 0, 'o'},
      {"runtime",             no_argument, 0,   0},
      {"outprefix",     required_argument, 0,   0},
      {"help",                no_argument, 0, 'h'},
      {"seed",          required_argument, 0, 's'},
      {"neventsperfile",required_argument, 0,   0},
      {0,                               0, 0,   0}
    };

    c = getopt_long(argc, argv, "o:n:s:h",
        long_options, &option_index);

    if( c == -1 ) break;

    switch (c) {
      case 'h':
        ShowHelp(argv[0]);
        exit(EXIT_SUCCESS);
        break;
      case 'o': SetOutputDirectory(std::string(optarg));break;
      case 'n': SetNumEvents(std::atoi(optarg)); break;
      case 's': SetRandomSeed(std::stoul(optarg)); break;
      case 0:
        switch (option_index) {
          case 0: SetFluxEnergyMin(std::atof(optarg)); break;
          case 1: SetFluxEnergyMax(std::atof(optarg)); break;
          case 2: SetFluxTimeMin(std::atof(optarg)); break;
          case 3: SetFluxTimeMax(std::atof(optarg)); break;
          case 4: SetRuntimeFactor(std::atof(optarg)); break;
          case 5: SetNumEvents(std::atoi(optarg)); break;
          case 6: SetOutputDirectory(std::string(optarg)); break;
          case 7: SetNormRuntime(true); break;
          case 8: SetOutputPrefix(std::string(optarg)); break;
          case 9: SetRandomSeed(std::stoul(optarg)); break;
          case 10: SetNumEventsPerFile(std::atoi(optarg)); break;
          default:
            ShowHelp(argv[0]);
            exit(EXIT_FAILURE);
            break;
        }
        break;
      default:
        ShowHelp(argv[0]);
        exit(EXIT_FAILURE);
        break;
    }
  }

  if (optind < argc ) {
    SetOutputDirectory(argv[optind]);
  }
}

bool SKSNSimUserConfiguration::CheckHealth () const {
  std::cout << "= Check configuration =" << std::endl;
  bool health = true;
  health &= CheckFluxEnergyMin();
  health &= CheckFluxEnergyMax();
  health &= CheckFluxTimeMin();
  health &= CheckFluxTimeMax();
  health &= CheckNormRuntime();
  health &= CheckNumEvents();
  health &= CheckRuntimeFactor();
  std::cout << "Config ==> good? " << health << std::endl;
  return health;
}

void SKSNSimUserConfiguration::Dump() const {
  std::cout << "= Dump configuration =" << std::endl;
  std::cout << "FluxEnergyMin = " << GetFluxEnergyMin() << std::endl;
  std::cout << "FluxEnergyMax = " << GetFluxEnergyMax() << std::endl;
  std::cout << "FluxTimeMin = " << GetFluxTimeMin() << std::endl;
  std::cout << "FluxTimeMax = " << GetFluxTimeMax() << std::endl;
  std::cout << "OutputDirecotry = " << GetOutputDirectory() << std::endl;
  std::cout << "OutputPrefix = " << GetOutputPrefix() << std::endl;
  std::cout << "NormRuntime = " << GetNormRuntime() << std::endl;
  std::cout << "NumEvents = " << GetNumEvents() << std::endl;
  std::cout << "NumEventsPerFile = " << GetNumEventsPerFile() << std::endl;
  std::cout << "RuntimeFactor = " << GetRuntimeFactor() << std::endl;
  std::cout << "RandomSeed = " << GetRandomSeed() << std::endl;
  std::cout << "====> Fine?  " << CheckHealth() << std::endl;

}

