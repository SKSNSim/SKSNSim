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

void SKSNSimUserConfiguration::ShowHelpDSNB(const char *argv0){
  std::cout << argv0
    << " [-c,--customflux {flux_filename}]"
    << " [--energy_min {energy_MeV}]"
    << " [--energy_max {energy_MeV}]"
    << " [--runtimefactor {float}]"
    << " [-n,--nevents {int}]"
    << " [-o,--outdir {directory}]"
    << " [--runtime]"
    << " [--runtime_begin {runnum}]"
    << " [--runtime_end {runnum}]"
    << " [--runtime_period {5/6}]"
    << " [--neventsperfile {int}]"
    << " [--outprefix {pref}]"
    << " [--flatposflux]"
    << " [-h,--help]"
    << " [-s,--seed {unsigned}]"
    << " [outputdirectory]"
    << std::endl << std::endl;

  std::cout << "Note: {} is essencial arguments, and [] is optional arguments" << std::endl << std::endl;
  std::cout << "Arguments: "  << std::endl
    << " -c,--customflux {flux_filename}: this option enforce to use specified flux file which should be formatted with \"energy(MeV) flux\" ( default = " << SKSNSimUserConfiguration::GetDefaultDSNBFluxModel() << " )" << std::endl
    << " --energy_min {energy_MeV}: lower energy limit to be generated in MeV ( default = " << SKSNSimUserConfiguration::GetDefaultFluxEnergyMin(SKSNSimUserConfiguration::MODEGENERATOR::kDSNB) << " MeV )" << std::endl
    << " --energy_max {energy_MeV}: uppwer energy limit to be generated in MeV ( default = " << SKSNSimUserConfiguration::GetDefaultFluxEnergyMax(SKSNSimUserConfiguration::MODEGENERATOR::kDSNB) << " MeV )" << std::endl
    << " --flatposflux: generate flat positron energy in range between --energy_min and --energy_max. " <<  std::endl
    << " --runtimefactor {float}: number of events per day for runtime normalization ( default = " << SKSNSimUserConfiguration::GetDefaultRuntimeNormFactor() << " evt/day )" << std::endl
    << " -n,--nevents {int}: number of events to be generated (exclusive with --runtime and --runtimefactor) ( default = " << SKSNSimUserConfiguration::GetDefaultNumEvents() << " )" << std::endl
    << " -o,--outdir {directory}: output directory. The generator fill events in the filename: {outdir}/{outprefix}_000000.root... ( default = " << SKSNSimUserConfiguration::GetDefaultOutputDirectory() << " )"  << std::endl
    << " --runtime: this option enebles runtime normalization. ( default = " << SKSNSimUserConfiguration::GetDefaultRuntimeNormalization() << " )"  << std::endl
    << " --runtime_begin {runnum}: begging run-number of runtime normalization ( default = " << SKSNSimUserConfiguration::GetDefaultRuntimeBegin() << " )"  << std::endl
    << " --runtime_end {runnum}: end run-nummber of runtime normalization ( default = " << SKSNSimUserConfiguration::GetDefaultRuntimeEnd() << " )" << std::endl
    << " --runtime_period {5/6}: SK period for runtime normalization: SK-IV = 3, SK-V = 4, SK-VI = 5 ( default = " << SKSNSimUserConfiguration::GetDefaultRuntimePeriod() << " , -1 means that period is specified by runnnumber. see --runtime_{begin/end}.)" << std::endl
    << " --neventsperfile {int}: number of events per one file ( default = " << SKSNSimUserConfiguration::GetDefaultNumEventsPerFile() << " )" << std::endl
    << " --outprefix {pref}: prefix of output file name ( default = " << SKSNSimUserConfiguration::GetDefaultOutputPrefix() << " )" << std::endl
    << " -h,--help: show this help" << std::endl
    << " -s,--seed {unsigned}: random seed ( default = " << SKSNSimUserConfiguration::GetDefaultRandomSeed() << " )" << std::endl
    << " {outputdirectory}: output direcotry. Same with -o,--outdir. This has higher priority."
    << std::endl
    << std::endl;
}

void SKSNSimUserConfiguration::ShowHelpSN(const char *argv0){

  std::cout << argv0
    << " [-h,--help]"
    << " [-o,--outdir outputdirectory]"
    << " [-m,--snmodel model_name]"
    << " [--nuosc 0(NONE)/1(NORMAL)/2(INVERTED)]"
    << " [-d,--distance distance_in_kpc]"
    << " [-g,--fillevent {0(no: just calculate expected num of evt)/1(yes: fill kinematics for detector sim.)}]"
    << " [--neventsperfile numberofevents]"
    << " [-s,--seed randomseed]"
    << " [--energy_min energy_MeV]"
    << " [--energy_max energy_MeV]"
    << " [--energy_nbins nbins]"
    << " [--runnum runnum]"
    << " [--subrunnum runnum]"
    << " [--time_min time]"
    << " [--time_max time]"
    << " [--time_nbins nbins]"
    << " [--outprefix prefix]"
    << " {outputdirectory}"
    << std::endl
    << std::endl;
  std::cout << "Note: {} is essencial arguments, and [] is optional arguments" << std::endl << std::endl;
  std::cout << "Usage: for old argument format" << std::endl
    << argv0
    << " {model_name}"
    << " {neutrino-oscillation}"
    << " {distance (normalized to 10kpc)}"
    << " {generate_event_or_not}"
    << " {output_directory}"
    << " {random_seed}"
    << std::endl << std::endl;
  std::cout << "Arguments: "  << std::endl
    << " -h,--help: show this help" << std::endl
    << " -o,--outdir {outputdirectory}: output directory (default = " << SKSNSimUserConfiguration::GetDefaultOutputDirectory() << " )" << std::endl
    << " -m,--snmodel {model_name}: name of SN flux model (default = " << SKSNSimUserConfiguration::GetDefaultSNModelName() << " ) " << std::endl
    << " --nuosc {int}: neutrino oscillation model: 0=NONE / 1=NORMAL / 2=INVERTED ( default = " << (int)SKSNSimUserConfiguration::GetDefaultNeutrinoOscType() << " )" << std::endl
    << " -d,--distance {distance_in_kpc}: distance from SN in unit of kpc ( default = " << SKSNSimUserConfiguration::GetDefaultSNDistanceKPC() << " kpc)" << std::endl
    << " -g,--fillevent [int]: if generate event kinematics for detector simulator: 0 = \"NO(just calculate expected num of evt) / 1 = YES (fill kinematics for detector sim.) (default = " << SKSNSimUserConfiguration::GetDefaultVectorGeneration() << " ). If you just specify \"-g\", turned ON" << std::endl
    << " --neventsperfile {numberofevents}: number of events in one file (default = " << SKSNSimUserConfiguration::GetDefaultNumEventsPerFile() << " )" << std::endl
    << " -s,--seed {randomseed}: random seed (default = " << SKSNSimUserConfiguration::GetDefaultRandomSeed() << " )" << std::endl
    << " --energy_min {energy_MeV}: energy lower limit in MeV (default = " << SKSNSimUserConfiguration::GetDefaultFluxEnergyMin(SKSNSimUserConfiguration::MODEGENERATOR::kSNBURST) << " MeV )" << std::endl
    << " --energy_max {energy_MeV}: energy upper limit in MeV (default = " << SKSNSimUserConfiguration::GetDefaultFluxEnergyMax(SKSNSimUserConfiguration::MODEGENERATOR::kSNBURST) << " MeV )" << std::endl
    << " --energy_nbins {nbins}: number of bins for energy (default = " << SKSNSimUserConfiguration::GetDefaultEnergyNBins() << " )" << std::endl
    << " --runnum {runnum}: run-number for MC (default =  none (999999))" << std::endl
    << " --subrunnum {runnum}: sub-run-number for MC (default = 0)" << std::endl
    << " --time_min {time_sec}: time lower limit in second (default = " << SKSNSimUserConfiguration::GetDefaultFluxTimeMin() << " sec )" << std::endl
    << " --time_max {time_sec}: time upper limit in second (default = " << SKSNSimUserConfiguration::GetDefaultFluxTimeMax() << " sec )" << std::endl
    << " --time_nbins {nbins}: number of bins for time (default = " << SKSNSimUserConfiguration::GetDefaultTimeNBins() << " )" << std::endl
    << " --outprefix {prefix}: prefix of output file name (default = " << SKSNSimUserConfiguration::GetDefaultOutputPrefix() << " )" << std::endl
    << std::endl;
  std::cout << "Arguments for old format"  << std::endl
    << " {model_name}: name of SN flux model" << std::endl
    << " {neutrino-oscillation}: neutrino oscillation model: 0=NONE / 1=NORMAL / 2=INVERTED" << std::endl
    << " {distance (normalized to 10kpc)}: distance from SN in unit of 10kpc" << std::endl
    << " {generate_event_or_not}: if generate event kinematics for detector simulator: 0 = \"NO(just calculate expected num of evt) / 1 = YES (fill kinematics for detector sim.)" << std::endl
    << " {output_directory}: output directory " << std::endl
    << " {random_seed}: random seed " << std::endl
    << std::endl;
}

void SKSNSimUserConfiguration::LoadFromArgsDSNB(int argc, char *argv[]){
  int c;
  while(1){
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;

    static struct option long_options[] = {
      {"customflux",    required_argument, 0, 'c'},
      {"energy_min",    required_argument, 0,   0},
      {"energy_max",    required_argument, 0,   0},
      {"runtimefactor", required_argument, 0,   0},
      {"nevents",       required_argument, 0, 'n'},
      {"outdir",        required_argument, 0, 'o'},
      {"runtime",             no_argument, 0,   0}, // 6
      {"runtime_begin", required_argument, 0,   0},
      {"runtime_end",   required_argument, 0,   0},
      {"runtime_period",required_argument, 0,   0},
      {"neventsperfile",required_argument, 0,   0},
      {"outprefix",     required_argument, 0,   0}, // 11 
      {"help",                no_argument, 0, 'h'},
      {"seed",          required_argument, 0, 's'},
      {"flatposflux",         no_argument, 0,   0},
      {0,                               0, 0,   0}
    };

    c = getopt_long(argc, argv, "c:n:o:hs:",
        long_options, &option_index);

    if( c == -1 ) break;


    switch (c) {
      case 'c': SetDSNBFluxModel(optarg); break;
      case 'h':
        ShowHelpDSNB(argv[0]);
        exit(EXIT_SUCCESS);
        break;
      case 'o': SetOutputDirectory(std::string(optarg));break;
      case 'n': SetNumEvents(std::atoi(optarg)); break;
      case 's': SetRandomSeed(std::stoul(optarg)); break;
      case 0:
        switch (option_index) {
          case 1: SetFluxEnergyMin(std::atof(optarg)); break;
          case 2: SetFluxEnergyMax(std::atof(optarg)); break;
          case 3: SetRuntimeFactor(std::atof(optarg)); break;
          case 6: SetNormRuntime(true); break;
          case 7: SetRuntimeRunBegin(std::atoi(optarg)); break;
          case 8: SetRuntimeRunEnd(std::atoi(optarg)); break;
          case 9: SetRuntimePeriod(std::atoi(optarg)); break;
          case 10: SetNumEventsPerFile(std::atoi(optarg)); break;
          case 11: SetOutputPrefix(optarg); break;
          case 14: SetDSNBFlatFlux(true); break;
          default:
            ShowHelpDSNB(argv[0]);
            exit(EXIT_FAILURE);
            break;
        }
        break;
      default:
        ShowHelpDSNB(argv[0]);
        exit(EXIT_FAILURE);
        break;
    }
  }

  if (optind < argc ) {
    SetOutputDirectory(argv[optind]);
  }
}

void SKSNSimUserConfiguration::LoadFromArgsSN(int argc, char *argv[]){
  int c;
  while(1){
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;

    static struct option long_options[] = {
      {"energy_min",    required_argument, 0,   0},
      {"energy_max",    required_argument, 0,   0},
      {"energy_nbins",  required_argument, 0,   0},
      {"time_min",      required_argument, 0,   0},
      {"time_max",      required_argument, 0,   0},
      {"time_nbins",    required_argument, 0,   0},
      {"snmodel",       required_argument, 0, 'm'}, // 6
      {"distance",      required_argument, 0, 'd'},
      {"fillevent",     optional_argument, 0, 'g'},
      {"nuosc",         required_argument, 0,   0},//9
      {"outdir",        required_argument, 0, 'o'},
      {"outprefix",     required_argument, 0,   0}, // 11
      {"help",                no_argument, 0, 'h'},
      {"seed",          required_argument, 0, 's'},
      {"neventsperfile",required_argument, 0,   0},
      {"runnum",        required_argument, 0,   0},
      {"subrunnum",     required_argument, 0,   0},
      {0,                               0, 0,   0}
    };

    c = getopt_long(argc, argv, "ho:m:d:g::s:",
        long_options, &option_index);

    if( c == -1 ) break;

    switch (c) {
      case 'h':
        ShowHelpSN(argv[0]);
        exit(EXIT_SUCCESS);
        break;
      case 'o': SetOutputDirectory(std::string(optarg));break;
      case 'm': SetSNBurstFluxModel(optarg); break;
      case 'd': SetSNDistanceKpc(std::stof(optarg)); break;
      case 'g': if( optarg == 0) SetVectorGeneration(true);
                  else SetVectorGeneration(std::atoi(optarg));
                  break;
      case 's': SetRandomSeed(std::stoul(optarg)); break;
      case 0:
        switch (option_index) {
          case 0: SetFluxEnergyMin(std::atof(optarg)); break;
          case 1: SetFluxEnergyMax(std::atof(optarg)); break;
          case 2: SetEnergyNBins(std::atoi(optarg)); break;
          case 3: SetFluxTimeMin(std::atof(optarg)); break;
          case 4: SetFluxTimeMax(std::atof(optarg)); break;
          case 5: SetTimeNBins(std::atoi(optarg)); break;
          case 9: SetNeutrinoOscType(std::atoi(optarg)); break;
          case 11: SetOutputPrefix(std::string(optarg)); break;
          case 14: SetNumEventsPerFile(std::atoi(optarg)); break;
          case 15: SetRunnum(std::atoi(optarg)); break;
          case 16: SetSubRunnum(std::atoi(optarg)); break;
          default:
            ShowHelpSN(argv[0]);
            exit(EXIT_FAILURE);
            break;
        }
        break;
      default:
        ShowHelpSN(argv[0]);
        exit(EXIT_FAILURE);
        break;
    }
  }

  /* check if old format or not */
  if (optind + 6 == argc ) {
    // OLD FORMAT
    SetSNBurstFluxModel(     argv[optind+0] );
    SetNeutrinoOscType( std::atoi(argv[optind+1]) );
    SetSNDistanceKpc( std::atof(argv[optind+2]) * 10. ); // in order to normalize 10 kpc, multiplying 10.
    SetVectorGeneration( std::atoi(argv[optind+3]) );
    SetOutputDirectory( argv[optind+4] );
    SetRandomSeed( std::atoi(argv[optind+5]) );
  }
  else if (optind < argc ) {
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
  std::cout << "FluxEnergyMin (MeV) = " << GetFluxEnergyMin() << std::endl;
  std::cout << "FluxEnergyMax (MeV) = " << GetFluxEnergyMax() << std::endl;
  std::cout << "EnergyNBins = " << GetEnergyNBins() << std::endl;
  std::cout << "FluxTimeMin (sec) = " << GetFluxTimeMin() << std::endl;
  std::cout << "FluxTimeMax (sec) = " << GetFluxTimeMax() << std::endl;
  std::cout << "TimeNBins = " << GetTimeNBins() << std::endl;
  std::cout << "OutputDirecotry = " << GetOutputDirectory() << std::endl;
  std::cout << "OutputPrefix = " << GetOutputPrefix() << std::endl;
  std::cout << "NumEventsPerFile = " << GetNumEventsPerFile() << std::endl;
  std::cout << "EventVectorGeneration = " << GetEventVectorGeneration() << std::endl;
  std::cout << "EventgenVolume = " << (int)GetEventgenVolume() << std::endl;
  std::cout << "NormRuntime = " << GetNormRuntime() << std::endl;
  std::cout << "NumEvents = " << GetNumEvents() << std::endl;
  std::cout << "RuntimeFactor = " << GetRuntimeFactor() << std::endl;
  std::cout << "RuntimeRunBegin = " << GetRuntimeRunBegin() << std::endl;
  std::cout << "RuntimeRunEnd = " << GetRuntimeRunEnd() << std::endl;
  std::cout << "RuntimePeriod = " << GetRuntimePeriod() << std::endl;
  std::cout << "Runnum = " << GetRunnum() << std::endl;
  std::cout << "SubRunnum = " << GetSubRunnum() << std::endl;
  std::cout << "SNDistance ( kpc ) = " << GetSNDistanceKpc() << std::endl;
  std::cout << "SNBurstFluxModel = " << GetSNBurstFluxModel() << std::endl;
  std::cout << "DSNBFluxModel = " << GetDSNBFluxModel() << std::endl;
  std::cout << "DSNBFlatFlux = " << GetDSNBFlatFlux() << std::endl;
  std::cout << "NuOscType = " << (int)GetNuOscType() << std::endl;
  std::cout << "RandomSeed = " << GetRandomSeed() << std::endl;
  std::cout << "====> Fine?  " << CheckHealth() << std::endl;

}

void SKSNSimUserConfiguration::Apply( SKSNSimVectorSNGenerator &gen ) const {
  gen.SetEnergyMin( GetFluxEnergyMin() );
  gen.SetEnergyMax( GetFluxEnergyMax() );
  gen.SetEnergyNBins( GetEnergyNBins() );
  gen.SetTimeMin( GetFluxTimeMin() );
  gen.SetTimeMax( GetFluxTimeMax() );
  gen.SetTimeNBins( GetTimeNBins() );
  gen.SetFlagFillEvent( GetEventVectorGeneration() );
  gen.SetGeneratorVolume( GetEventgenVolume() );
  gen.SetSNDistanceKpc( GetSNDistanceKpc() );
  gen.SetGeneratorNuOscType( GetNuOscType() );
  gen.SetRUNNUM( GetRunnum() );
  gen.SetSubRUNNUM( GetSubRunnum() );
  gen.SetRandomSeed( GetRandomSeed() );
  std::cout << "getTiemNBins= " << GetTimeNBins() << std::endl;
}

void SKSNSimUserConfiguration::Apply( SKSNSimVectorGenerator &gen ) const {
  gen.SetEnergyMin( GetFluxEnergyMin() );
  gen.SetEnergyMax( GetFluxEnergyMax() );
  gen.SetGeneratorVolume( GetEventgenVolume() );
  gen.SetNormRuntime( GetNormRuntime() );
  gen.SetRuntimeFactor( GetRuntimeFactor() );
  gen.SetRuntimeBegin( GetRuntimeRunBegin() );
  gen.SetRuntimeEnd( GetRuntimeRunEnd() );
  gen.SetRuntimePeriod( GetRuntimePeriod() );
  gen.SetFlatPositronFlux( GetDSNBFlatFlux() );
  gen.SetRUNNUM( GetRunnum() );
  gen.SetSubRUNNUM( GetSubRunnum() );
}

SKSNSimUserConfiguration::MODERUNTIME SKSNSimUserConfiguration::CheckMODERuntime() const {
  if( GetNormRuntime() ){
    if( GetRuntimePeriod() == -1 )
      return SKSNSimUserConfiguration::MODERUNTIME::kRUNTIMERUNNUM;
    else
      return SKSNSimUserConfiguration::MODERUNTIME::kRUNTIMEPERIOD;
  }

  return SKSNSimUserConfiguration::MODERUNTIME::kEVNUM;

}
