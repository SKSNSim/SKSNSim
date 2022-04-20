
#include "VectGenIO.hh"

void ShowUsage(const int argc, char **argv){
  std::cerr<<" Usage: " << argv[0] << " [options] {} {} {}" << std::endl;
  std::cerr<<" Argument for "<< argv[0] << std::endl;
  std::cerr<<" 1st: Number of event" << std::endl;
  std::cerr<<" 2nd: Output directory, default: ./data/"<< std::endl;
  std::cerr<<" 3rd: Random seed"<< std::endl;
  std::cerr<<" Options: -f\n"
    << "    -f {flux_file_name}: file name of flux data"
    << std::endl;
  std::cerr<<" return -1"<<std::endl;
}

int main( int argc, char ** argv )
{

	std::cout<<" Start Generating Event "<< argc <<std::endl;

	if(argc < 3){
    ShowUsage(argc, argv);
	  exit(1);
	}

	/*-----Number of Generated Event-----*/
	int NumEv = atoi(argv[1]);

	/*-----Output directory-----*/
	std::string OutDirIO("./data/");
	OutDirIO = std::string(argv[2]);

	/*-----random number initialization-----*/
	uint seedIO = 0;
	seedIO = atoi(argv[3]);

  /*---- Flux data ----------------*/
  std::string FluxFile("../expect/horiuchi/8MeV_Nominal.dat");
	//FluxFile = std::string(argv[4]);


	/*-----Geneartion-----*/

	VectGenIO *io = new VectGenIO(OutDirIO, seedIO, FluxFile);

	io->DoProcess(NumEv);

}
