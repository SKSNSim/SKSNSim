
#include "VectGenIO.hh"

int main( int argc, char ** argv )
{

	std::cout<<" Start Generating Event "<< argc <<std::endl;

	if(argc < 3){
	  std::cerr<<" Argument for "<< argv[0] << std::endl;
	  std::cerr<<" 1st: Number of event" << std::endl;
	  std::cerr<<" 2nd: Output directory, default: ./data/"<< std::endl;
	  std::cerr<<" 3rd: Random seed"<< std::endl;
	  std::cerr<<" return -1"<<std::endl;
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

	/*-----Geneartion-----*/
	VectGenIO *io = new VectGenIO(OutDirIO, seedIO);

	io->DoProcess(NumEv);

}
