
#include "VectGenIO.hh"

int main( int argc, char ** argv )
{
	std::cout<<" Start Generating Event "<<std::endl;

	if(argc < 2){
		std::cerr<<"Usage: "<< argv[0]
			<<" model name (e.g. nakazato/intp2001.data) neutrino-oscillation (0:no, 1:normal, 2:inverted)"<<std::endl;
		std::cerr<<"return -1"<<std::endl;
		exit(1);
	}

	/*-----set model name-----*/
	std::string ModelName = std::string(argv[1]);

	/*-----neutrino oscillation-----*/
	int nuosc = atoi(argv[2]);

	/*-----Geneartion-----*/
	VectGenIO *io = new VectGenIO(ModelName, nuosc);
	io->DoProcess();

}
