#include "SKSNSimVectorGenerator.hh"
#include <memory>
#include <TRandom3.h>
#include "SKSNSimTools.hh"
#include "SKSNSimFileIO.hh"
#include "SKSNSimUserConfiguration.hh"

int main( int argc, char ** argv )
{

	//std::cout<<" Start Generating Event "<< argc <<std::endl;

	//if(argc < 5){
	//  std::cerr<<" Argument for "<< argv[0] << std::endl;
	//  std::cerr<<" 1st: model name (e.g. nakazato/intp2001.data)" << std::endl;
	//  std::cerr<<" 2nd: neutrino-oscillation (0:no, 1:normal, 2:inverted)"<< std::endl;
	//  std::cerr<<" 3rd: Distance (normalized to 10kpc)"<< std::endl;
	//  std::cerr<<" 4th: Generate event or not (0 : No just calculate expected, 1: Yes for detector simulation)"<< std::endl;
	//  std::cerr<<" 5th: Output directory, default: ./data/(SN model)"<< std::endl;
	//  std::cerr<<" 6th: Random seed"<< std::endl;
	//  std::cerr<<" return -1"<<std::endl;
	//  exit(1);
	//}

	///*-----set model name-----*/
	//std::string ModelName = std::string(argv[1]);

	///*-----neutrino oscillation-----*/
	//int nuosc = atoi(argv[2]);

	///*-----Distance-----*/
	//double distIO = atof(argv[3]); // normalized to 10kpc

	///*-----Generate events or not-----*/
	//int flagIO = atoi(argv[4]);

	///*-----Output directory-----*/
	//std::string OutDirIO("./data");
	//if(argc > 5) OutDirIO = std::string(argv[5]);

	///*-----random number initialization-----*/
	//uint seedIO = 0;
	//if(argc > 6) seedIO = atoi(argv[6]);
  //
  //
  auto config = std::make_unique<SKSNSimUserConfiguration>();
  config->LoadFromArgsSN(argc, argv);
  config->Dump();


  std::unique_ptr<TRandom3> rng = std::make_unique<TRandom3>(42);
  std::unique_ptr<SKSNSimSNFluxNakazatoFormat> flux (new SKSNSimSNFluxNakazatoFormat());
  flux->SetModel(config->GetSNBurstFluxModel());
  

	/*-----Geneartion-----*/
  std::unique_ptr<SKSNSimVectorSNGenerator> generator = std::make_unique<SKSNSimVectorSNGenerator>();
  generator->SetRandomGenerator(std::move(rng));
  generator->AddFluxModel(std::move(flux));
  config->Apply(*generator);
  auto buffer = generator->GenerateEvents();
  SKSNSimTools::DumpDebugMessage(Form("Successed GenerateEvents -> %d events", (int)buffer.size()));

  auto flist = GenerateOutputFileList(*config);
  for(auto it = flist.begin(); it != flist.end(); it++){
    auto vectio = std::make_unique<SKSNSimFileOutTFile>();
    vectio->Open(it->GetFileName(), true);
    vectio->Write(buffer);
    vectio->Close();
  }

  return EXIT_SUCCESS;
}
