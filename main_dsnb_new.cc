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
#include "SKSNSimCrosssection.hh"


int main(int argc, char **argv){

  auto vectio = std::make_unique<SKSNSimFileOutTFile>("test.root");
  auto vectgen = std::make_unique<SKSNSimVectorGenerator>();
  auto fluxhoriuch = std::make_unique<SKSNSimFluxDSNBHoriuchi>();
  fluxhoriuch->DumpFlux();
  auto xsec = std::make_unique<SKSNSimXSecIBDSV>();
  vectgen->SetRandomSeed(42);
  vectgen->AddFluxModel((SKSNSimFluxModel*)fluxhoriuch.release());
  vectgen->AddXSecModel((SKSNSimCrosssectionModel*)xsec.release());
  vectgen->SetEnergyMin(3.8);
  vectgen->SetEnergyMax(50.0);

  for(int i = 0; i < 10000; i++){
    vectio->Write(vectgen->GenerateEventIBD());
  }
  vectio->Close();



  return EXIT_SUCCESS;
}
