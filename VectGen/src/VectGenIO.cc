
#include "VectGenSetBin.hh"
#include "VectGenIO.hh"

VectGenIO::VectGenIO(std::string ModelName, int nuosc)
//  : Generator()
{

	/*-----set binning-----*/
	VectGenSetBinValues();

	/*-----set neutrino oscillation parameter -----*/
	switch(nuosc){
	case 0:
		oscnue1 = 1.;
		oscnue2 = 0.;
		oscneb1 = 1.;
		oscneb2 = 0.;
		oscnux1 = 2.; // nu_mu and nu_tau
		oscnux2 = 0.;
		oscnxb1 = 2.; // nu_mu_bar and nu_tau_bar
		oscnxb2 = 0.;
		break;

	case 1:
		oscnue1 = 0.;
		oscnue2 = 1.;
		oscneb1 = cos2th12;
		oscneb2 = sin2th12;
		oscnux1 = 1.;
		oscnux2 = 1.;
		oscnxb1 = 1.+ cos2th12;
		oscnxb2 = sin2th12;
		break;

	case 2:
		oscnue1 = sin2th12;
		oscnue2 = cos2th12;
		oscneb1 = 0.;
		oscneb2 = 1.;
		oscnux1 = 1.+ sin2th12;
		oscnux2 = cos2th12;
		oscnxb1 = 1.;
		oscnxb2 = 1.;
		break;

	default :
		std::cout << "3rd argument should be 0 or 1 or 2, but it is " << nuosc << std::endl;
		exit(1);
	}

	totNuebarp=0.;
	totNueElastic=0., totNuebarElastic=0., totNuxElastic=0., totNuxbarElastic=0.;
	totNueO=0., totNuebarO=0.;
	totNcNup=0., totNcNun=0., totNcNubarp=0., totNcNubarn=0.;

	//Class call of neutrino flux table
	std::stringstream ssname1;
	ssname1 << "/home/sklowe/supernova/data/nakazato/" << ModelName;
	std::string FileIn = ssname1.str();
	nuflux = new VectGenSnNakazato(FileIn.c_str());
	//nuflux = new SnWilson();

	//Class call of cross-section calculation
	nucrs = new VectGenNuCrosssection(); //changed from VB to Vissani
	nucrs->ReadCsNuElastic();


	// set root file name
	std::stringstream ssname2;
	ssname2 << "./data/" << ModelName << ".root";
	FileRoot = ssname2.str();
	std::cout << "root file name " << FileRoot << std::endl;


	// set text file name
	std::stringstream ssname3;
	ssname3 << "./data/"<< ModelName << ".dat";
	FileText = ssname3.str();
	std::cout << "text file name " << FileText << std::endl;

	// set histogram
	evrate = new TH1D("evrate","", 100, 0., 100.);

}

void VectGenIO::DoProcess()
{
  Process();

  TFile* fout = new TFile(FileRoot.c_str(), "RECREATE");
  evrate->Write();
  fout->Close();

  std::cout << "each event number is "<< totNuebarp << " " << totNueElastic << " " << totNuebarElastic << " " << totNuxElastic << " " << totNuxbarElastic << std::endl;

  std::ofstream ofs(FileText.c_str());
  ofs << totNuebarp << " " << totNueElastic << " " << totNuebarElastic << " " << totNuxElastic << " " << totNuxbarElastic << std::endl;

}
