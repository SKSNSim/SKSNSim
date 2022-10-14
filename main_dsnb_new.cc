/******************************
 * File: main_dsnb_new.cc
 * Description:
 * new vector generator for DSNB simulation
 *******************************/

#include <string>
#include <memory>
#include "SKSNSimVectorGenerator.hh"
#include "SKSNSimFileIO.hh"
#include "SKSNSimFlux.hh"
#include "SKSNSimTools.hh"
#include "SKSNSimCrosssection.hh"
#include "SKSNSimUserConfiguration.hh"


int main(int argc, char **argv){

  auto config = std::make_unique<SKSNSimUserConfiguration>();
  config->LoadFromArgsDSNB(argc, argv);
  config->Dump();
  config->CheckHealth();

  auto vectgen = std::make_unique<SKSNSimVectorGenerator>();
  config->Apply(*vectgen);
  auto fluxhoriuch = std::make_unique<SKSNSimFluxDSNBHoriuchi>();
  //fluxhoriuch->DumpFlux();
  auto xsec = std::make_unique<SKSNSimXSecIBDSV>();

  vectgen->SetRandomGenerator(config->GetRandomGenerator());
  vectgen->AddFluxModel((SKSNSimFluxModel*)fluxhoriuch.release());
  vectgen->AddXSecModel((SKSNSimCrosssectionModel*)xsec.release());

  if( config->CheckMODERuntime() == SKSNSimUserConfiguration::MODERUNTIME::kEVNUM ){
    auto flist = GenerateOutputFileList(*config);
    for(auto it = flist.begin(); it != flist.end(); it++){
      auto vectio = std::make_unique<SKSNSimFileOutTFile>();
      vectio->Open(it->GetFileName());
      size_t pos = std::distance(flist.begin(), it);
      for(int i = 0; i < it->GetNumEvents(); i++)
        vectio->Write(vectgen->GenerateEventIBD());
      vectio->Close();
    }
  } else {
    std::vector<std::tuple<int,double>> livetimes;
    if( config->CheckMODERuntime() == SKSNSimUserConfiguration::MODERUNTIME::kRUNTIMERUNNUM ){
      livetimes = SKSNSimLiveTime::LoadLiveTime( config->GetRuntimeRunBegin(), config->GetRuntimeRunEnd());
    } else if( config->CheckMODERuntime() == SKSNSimUserConfiguration::MODERUNTIME::kRUNTIMEPERIOD ){
      livetimes = SKSNSimLiveTime::LoadLiveTime( (SKSNSIMENUM::SKPERIOD) config->GetRuntimePeriod());
    }

  }

  return EXIT_SUCCESS;
}
