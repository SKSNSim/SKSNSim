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
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>
#include <mcinfo.h>
#include <SKSNSimCrosssection.hh>
#include <SKSNSimFlux.hh>


void showUsage( const char *argv0){
  std::cout << "Usage: " << argv0 << " [options] {output.root} {input.root}" << std::endl;
  std::cout << std::endl;
  std::cout << "Options: -h,--help" << std::endl
    << "-h,--help: show this help" << std::endl
    << std::endl;
}


int main( int argc, char **argv){

  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"help",     no_argument,      0, 'h'},
      {0,         0,                 0,  0 }
    };

    const int c = getopt_long(argc, argv, "h",
        long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
      case 0:
        printf("option %s", long_options[option_index].name);
        if (optarg)
          printf(" with arg %s", optarg);
        printf("\n");
        break;

      case '0':
      case '1':
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
  double flux_horiuchi_8mev;

  tr->SetBranchAddress("MC", &MC);
  auto newBrIBDXsec = tr->Branch("ibd_xsec", &ibd_xsec, "ibd_xsec/D");
  auto newBrFluxHoriuchi8MeV = tr->Branch("flux_horiuchi_8mev", &flux_horiuchi_8mev, "flux_horiuchi_8mev/D");

  const int ntotal = tr->GetEntries();


  // For xsec
  std::unique_ptr<SKSNSimXSecIBDRVV> xsec_rvv = std::make_unique<SKSNSimXSecIBDRVV>();

  // For flux
  std::unique_ptr<SKSNSimFluxDSNBHoriuchi> fluxcal_horiuchi_8mev = std::make_unique<SKSNSimFluxDSNBHoriuchi>();

  for( int e = 0; e < ntotal; e++){
    tr->GetEntry(e);

    const double nue_energy = MC->energy[0];
    const TVector3 nue_p ( MC->pvc[0][0], MC->pvc[0][1], MC->pvc[0][2]);
    const TVector3 pos_p ( MC->pvc[2][0], MC->pvc[2][1], MC->pvc[2][2]);

    ibd_xsec = xsec_rvv->GetDiffCrosssection( nue_energy, TMath::Cos(nue_p.Angle(pos_p))).first;
    flux_horiuchi_8mev = fluxcal_horiuchi_8mev->GetFlux(nue_energy, 0, SKSNSimFluxModel::FLUXNUTYPE::FLUXNUEB);


    newBrIBDXsec->Fill();
    newBrFluxHoriuchi8MeV->Fill();

  }

  tr->Write(nullptr, TObject::kWriteDelete);
  ofile->Save();
  ofile->Close();
  ifile->Close();

  return EXIT_SUCCESS;
}
