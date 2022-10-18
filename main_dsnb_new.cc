/******************************
 * File: main_dsnb_new.cc
 * Description:
 * new vector generator for DSNB simulation
 *******************************/

#include <string>
#include <memory>
#include <sstream>
#include <iomanip>
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

  auto flist = GenerateOutputFileList(*config);
  for(auto it = flist.begin(); it != flist.end(); it++){
    auto vectio = std::make_unique<SKSNSimFileOutTFile>(it->GetFileName());
    vectgen->SetRUNNUM( it->GetRun() );
    vectgen->SetSubRUNNUM( it->GetSubrun() );
    vectio->Write(vectgen->GenerateEvents(it->GetNumEvents()));
    vectio->Close();
  }

  return EXIT_SUCCESS;
}
