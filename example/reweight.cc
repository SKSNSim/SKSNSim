/*==================================================
 * reweight.cc
 *
 * Description:
 * Reweight the events according to generated IBD topology
 * and IBD cross section implemented in the SKSNSim
 *
 * Date: Fri Nov 10 16:46:59 JST 2023
 * Author: Shota Izumiyaam (izumiyama@hep.phys.titech.ac.jp)
 *
 * ==============================================*/

#include <iostream>
#include <memory>
#include <getopt.h>
#include <cstdlib>
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>
#include <mcinfo.h>
#include <SKSNSimCrosssection.hh>
#include <SKSNSimFlux.hh>

using namespace boost::filesystem;
using std::cout;

namespace PID {
  constexpr int NUE = 12;
  constexpr int ELECTRON = 11;
}


void showUsage( const char *argv0){
  std::cout << "Usage: " << argv0 << " [options] {output.root} {input.root}" << std::endl;
  std::cout << std::endl;
  std::cout << "Options: -h,--help,-d,--fluxdir, --fluxpositrondir" << std::endl
    << "-h,--help: show this help" << std::endl
    << "-d {directory}, --fluxdir {directory}: search for flux file (.txt) in directory {directory}" << std::endl
    << "--fluxpositrondir {directory}: search for flux file for positron spectrum (.txt) in directory {directory} for Li9" << std::endl
    << std::endl;
}

typedef std::pair<std::string, std::string> MODEL;
std::vector<MODEL> findModelFile(std::string directory) {
  std::vector<MODEL> buffer;

  if( ! is_directory(directory)) return buffer;
  for( directory_entry&x : directory_iterator(directory) ){
    if( ! is_regular_file(x) ) continue;
    if( x.path().extension() == ".txt" ) {
      std::string mname = x.path().stem().string();
      buffer.push_back( std::make_pair( mname, x.path().string()));
    }
  }
  return buffer;
}

class FLUXGENERATOR {
  public:
    std::string name;
    TBranch *br;
    Double_t flux;
    std::unique_ptr<SKSNSimDSNBFluxCustom> generator;
    FLUXGENERATOR(): br(nullptr) {
      //generator = std::make_unique<SKSNSimDSNBFluxCustom>();
    }
    FLUXGENERATOR(std::string n, std::string fname, TTree &tr): name(n) {
      TString brname = Form("flux_%s", name.c_str());
      br = tr.Branch( brname, &flux, brname+"/D");
      generator = std::make_unique<SKSNSimDSNBFluxCustom>(fname, " ");
      //cout << "Cnstr' " << name << " " << generator->GetFlux(15.0)<<std::endl;
    }
    void fill( double energy){
      flux = generator->GetFlux( energy);
      //cout << name << " " << flux << std::endl;
      br->Fill();
    }
};


int main( int argc, char **argv){

  std::vector<MODEL> flux_models;
  std::vector<MODEL> flux_models_positron;

  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"help",                  no_argument,      0, 'h'},
      {"fluxdir",         required_argument,      0, 'd'},
      {"fluxpositrondir", required_argument,      0,  0 },
      {0,                 0,                      0,  0 }
    };

    const int c = getopt_long(argc, argv, "d:h",
        long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
      case 0:
        switch( option_index) {
          case 2: // fluxpositrondir
            flux_models_positron = findModelFile( optarg);
            cout << "Loaded flux files (positron): " << flux_models_positron.size() << " files -> ";
            for( auto it : flux_models_positron ) cout << it.first << " ";
            cout << std::endl;
            break;
          default:
            break;
        }
        break;
      case 'd':
        flux_models = findModelFile( optarg);
        cout << "Loaded flux files: " << flux_models.size() << " files -> ";
        for( auto it : flux_models ) cout << it.first << " ";
        cout << std::endl;
        break;
      case 'h':
        showUsage(argv[0]);
        break;
      default:
        printf("?? getopt returned character code 0%o ??\n", c);
    }
  }

  if ( argc - optind != 2 ) {
    showUsage(argv[0]);
    exit(EXIT_FAILURE);
  }

  const std::string ofname = argv[optind++];
  const std::string ifname = argv[optind++];


  std::unique_ptr<TFile> ifile = std::make_unique<TFile>(ifname.c_str(), "READ");
  std::unique_ptr<TFile> ofile = std::make_unique<TFile>(ofname.c_str(), "RECREATE");
  ofile->cd();

  TTree *tr = ((TTree*) ifile->Get("data"))->CloneTree();

  MCInfo *MC = nullptr;
  double ibd_xsec;

  tr->SetBranchAddress("MC", &MC);
  auto newBrIBDXsec = tr->Branch("ibd_xsec", &ibd_xsec, "ibd_xsec/D");

  const int ntotal = tr->GetEntries();


  // For xsec
  std::unique_ptr<SKSNSimXSecIBDRVV> xsec_rvv = std::make_unique<SKSNSimXSecIBDRVV>();

  // For flux
  std::vector<std::unique_ptr<FLUXGENERATOR>> flux_container;
  for( auto it : flux_models)
    flux_container.push_back( std::make_unique<FLUXGENERATOR>( it.first, it.second, *tr));

  std::vector<std::unique_ptr<FLUXGENERATOR>> flux_container_positron;
  for( auto it : flux_models_positron)
    flux_container_positron.push_back( std::make_unique<FLUXGENERATOR>( it.first, it.second, *tr));



  for( int e = 0; e < ntotal; e++){
    tr->GetEntry(e);

    int positronID = -1;
    int nuebarID = -1;
    for( int i = 0; i < MC->nvc; i++){
      if( MC->ipvc[i] == -PID::NUE) nuebarID = i;
      else if( MC->ipvc[i] == -PID::ELECTRON) positronID = i;
    }
    const double nue_energy = nuebarID >= 0 ? MC->energy[nuebarID]: 0.0;
    const double positron_energy = positronID >= 0 ? MC->energy[positronID]: 0.0;
    const TVector3 nue_p = nuebarID >= 0 ? TVector3( MC->pvc[nuebarID][0], MC->pvc[nuebarID][1], MC->pvc[nuebarID][2]) : TVector3();
    const TVector3 pos_p = positronID >= 0 ? TVector3( MC->pvc[positronID][0], MC->pvc[positronID][1], MC->pvc[positronID][2]): TVector3();

    ibd_xsec = xsec_rvv->GetDiffCrosssection( nue_energy, TMath::Cos(nue_p.Angle(pos_p))).first;

    newBrIBDXsec->Fill();
    for( auto &it : flux_container) 
      it->fill( nue_energy);
    for( auto &it : flux_container_positron) 
      it->fill( positron_energy);
  }

  tr->Write(nullptr, TObject::kWriteDelete);
  ofile->Save();
  ofile->Close();
  ifile->Close();

  return EXIT_SUCCESS;
}
