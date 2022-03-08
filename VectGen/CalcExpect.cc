
#include "VectGenIO.hh"

int main( int argc, char ** argv )
{
  std::cout<<" Start Generating Event "<< argc <<std::endl;

	if(argc < 5){
	  std::cerr<<" Argument for "<< argv[0] << std::endl;
	  std::cerr<<" 1st: model name (e.g. nakazato/intp2001.data)" << std::endl;
	  std::cerr<<" 2nd: neutrino-oscillation (0:no, 1:normal, 2:inverted)"<< std::endl;
	  std::cerr<<" 3rd: Distance (normalized to 10kpc)"<< std::endl;
	  std::cerr<<" 4th: make event or not (0 : no, just calculate expected, 1: yes, for detector simulation)"<< std::endl;
	  std::cerr<<" return -1"<<std::endl;
	  exit(1);
	}

	/*-----set model name-----*/
	std::string ModelName = std::string(argv[1]);

	/*-----neutrino oscillation-----*/
	int nuosc = atoi(argv[2]);

	/*-----Distance-----*/
	double Distance = atof(argv[3]); // normalized to 10kpc

	/*-----making events or not-----*/
	int flag_event = atoi(argv[4]);

	/*-----Geneartion-----*/
	VectGenIO *io = new VectGenIO(ModelName, nuosc, flag_event);

	io->DoProcess(flag_event);

}
