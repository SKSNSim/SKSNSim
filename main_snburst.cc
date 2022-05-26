#include "VectGenIO.hh"

int main( int argc, char ** argv )
{

	std::cout<<" Start Generating Event "<< argc <<std::endl;

	if(argc < 5){
	  std::cerr<<" Argument for "<< argv[0] << std::endl;
	  std::cerr<<" 1st: model name (e.g. nakazato/intp2001.data)" << std::endl;
	  std::cerr<<" 2nd: neutrino-oscillation (0:no, 1:normal, 2:inverted)"<< std::endl;
	  std::cerr<<" 3rd: Distance (normalized to 10kpc)"<< std::endl;
	  std::cerr<<" 4th: Generate event or not (0 : No just calculate expected, 1: Yes for detector simulation)"<< std::endl;
	  std::cerr<<" 5th: Output directory, default: ./data/(SN model)"<< std::endl;
	  std::cerr<<" 6th: Random seed"<< std::endl;
	  std::cerr<<" return -1"<<std::endl;
	  exit(1);
	}

	/*-----set model name-----*/
	std::string ModelName = std::string(argv[1]);

	/*-----neutrino oscillation-----*/
	int nuosc = atoi(argv[2]);

	/*-----Distance-----*/
	double distIO = atof(argv[3]); // normalized to 10kpc

	/*-----Generate events or not-----*/
	int flagIO = atoi(argv[4]);

	/*-----Output directory-----*/
	std::string OutDirIO("./data/");
	if(argc > 5) OutDirIO = std::string(argv[5]);

	/*-----random number initialization-----*/
	uint seedIO = 0;
	if(argc > 6) seedIO = atoi(argv[6]);

	/*-----Geneartion-----*/
	VectGenIO *io = new VectGenIO(ModelName, nuosc, distIO, flagIO, OutDirIO, seedIO);

	io->DoProcess();

}
