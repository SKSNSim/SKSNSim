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

  /* Load configuration from arguments */
  auto config = std::make_unique<SKSNSimUserConfiguration>();
  config->LoadFromArgsDSNB(argc, argv);
  config->Dump();
  config->CheckHealth();

  /* Configure generator class */
  auto vectgen = std::make_unique<SKSNSimVectorGenerator>();
  config->Apply(*vectgen);
  vectgen->SetRandomGenerator(config->GetRandomGenerator());

  /* Setup of flux model and xsec model */
  auto flux = std::make_unique<SKSNSimDSNBFluxCustom>( config->GetDSNBFluxModel() );
  auto xsec = std::make_unique<SKSNSimXSecIBDSV>();
  vectgen->AddFluxModel((SKSNSimFluxModel*)flux.release());
  vectgen->AddXSecModel((SKSNSimCrosssectionModel*)xsec.release());

  /*  Tempolary variables to define integration of dN/dE spectrum */
  int num_random_throw = 0;
  int num_total_event = 0;
  double max_weight = 0;
  auto calcToralRandomThrow = [](const std::vector<SKSNSimSNEventVector> &v){
    int n = 0;
    for(auto it = v.begin(); it != v.end(); it++) n += it->GetNRandomThrow();
    return n;
  };


  /* Generate number of events and output file-name etc. */
  auto flist = GenerateOutputFileList(*config);

  for(auto it = flist.begin(); it != flist.end(); it++){
    /* Open file IO and then generate events from file configuration */
    auto vectio = std::make_unique<SKSNSimFileOutTFile>(it->GetFileName());
    vectgen->SetRUNNUM( it->GetRun() );
    vectgen->SetSubRUNNUM( it->GetSubrun() );
    auto evt_buffer = vectgen->GenerateEvents(it->GetNumEvents());

    /*  Calculate event weight in order to define integration of dN/dE spectrum */
    num_random_throw += calcToralRandomThrow(evt_buffer);
    num_total_event += evt_buffer.size();
    if( evt_buffer.size() > 0 && max_weight < evt_buffer.at(0).GetWeightMaxProb() ) max_weight = evt_buffer.at(0).GetWeightMaxProb();

    /* Output generated vectors and close file IO */
    vectio->Write(evt_buffer);
    vectio->Close();
  }

  /* Calculation of integrateion of dN/dE spectrum */
  if( num_random_throw > 0 ){
  std::cout << "=============================" << std::endl
    <<  "Finished event generation: integration results: " << std::endl
    << "(Total Events) / (Total Random Throw)  = " << num_total_event << " / " << num_random_throw << " = " << (double)num_total_event/(double)num_random_throw << std::endl
    << "Weight of max-probability in hit-and-miss method = " << max_weight << std::endl
    << "(Integration of dN/dE spectrum (flux x xsec)) / ( total number of free-proton ) = " << max_weight * (double) num_total_event / (double)num_random_throw << std::endl
    << "============================="   << std::endl;
  } else {
    std::cout <<  "Finished event generation: total number of random throw is zero or negative ( " << num_random_throw << " )" << std::endl
      << "Ignore this message if you do NOT use flat-flux mode." << std::endl;
  }

  return EXIT_SUCCESS;
}
