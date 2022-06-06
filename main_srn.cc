
#include "VectGenIO.hh"
#include "VectGenSetBin.hh"
#include <getopt.h>

constexpr double DEFAULT_NUENE_MIN = 10.0; // number from VectGenSetBin.hh::nuEneMin/Max
constexpr double DEFAULT_NUENE_MAX = 50.0;

void ShowUsage(char *arg0){
  std::cerr<< arg0 << " [options] {num_eve} {out_dir} {seed}" << std::endl;
  std::cerr<< "Usage: vector generetaor for DSBN MC" << std::endl;
  std::cerr<<" Argument for "<< arg0 << std::endl;
  std::cerr<<" 1st: Number of event" << std::endl;
  std::cerr<<" 2nd: Output directory, default: ./data/"<< std::endl;
  std::cerr<<" 3rd: Random seed"<< std::endl;
  std::cerr<<" Options: -f, -h, --flux_file, --nu_ene_{min,max} \n"
    << "\t-h, --help: show this help" << std::endl
    << "\t-f {flux_data_file}, --flux_file {flux_data_file}: specify flux data file to be used" << std::endl
    << "\t--nu_ene_min {energy_in_MeV}: minimum total energy [MeV] to be generated (default = " << DEFAULT_NUENE_MIN << " MeV)" << std::endl
    << "\t--nu_ene_max {energy_in_MeV}: maximum total energy [MeV] to be generated (default = " << DEFAULT_NUENE_MAX << " MeV)" << std::endl;
  std::cerr<<" return -1"<<std::endl;
}

int main( int argc, char ** argv )
{

	std::cout<<" Start Generating Event "<< argc <<std::endl;

  //=====================
  // argument controll
  int c, digit_optind = 0;
  double nuEne_min = -1.0;
  double nuEne_max = -1.0;
  std::string FluxFile ("../expect/horiuchi/8MeV_Nominal.dat");
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    unsigned long seed_read = 0;
    enum OPTIONS {kNuEneMin = 0, kNuEneMax, kFluxFile, kHelp};
    static struct option long_options[] = {
      {"nu_ene_min",  required_argument, 0, 0},
      {"nu_ene_max",  required_argument, 0, 0},
      {"flux_file",   required_argument, 0, 0},
      {"help",              no_argument, 0, 0},
      {0,                             0, 0, 0}
    };

    c = getopt_long( argc, argv, "hf:v",
        long_options, &option_index);
    if(c==-1)
      break;

    switch (c) {
      case 0:
        switch(option_index) {
          case kNuEneMin:
            nuEne_min = strtod(optarg, NULL);
            std::cerr << "GETOPT: minumum neutrino energy <- set to = " << nuEne_min << std::endl;
            break;
          case kNuEneMax:
            nuEne_max = strtod(optarg, NULL);
            std::cerr << "GETOPT: maximum neutrino energy <- set to = " << nuEne_max << std::endl;
            break;
          case kFluxFile:
            FluxFile = std::string(optarg);
            std::cerr << "GETOPT: flux data file = " << FluxFile << std::endl;
            break;
          case kHelp:
            ShowUsage(argv[0]);
            return EXIT_SUCCESS;
            break;
          default:
            break;
        }
        break;
      case 'f':
        FluxFile = std::string(optarg);
        std::cerr << "GETOPT: flux data file = " << FluxFile << std::endl;
        break;
      case 'h':
        ShowUsage(argv[0]);
        return EXIT_SUCCESS;
        break;
      default:
        std::cerr << "GETOPT: strange argument: code = " << c << std::endl;
    }
  }

  if( optind+3 > argc){
    ShowUsage(argv[0]);
    return EXIT_FAILURE;
  }

	/*-----Number of Generated Event-----*/
	int NumEv = atoi(argv[optind]);

	/*-----Output directory-----*/
	std::string OutDirIO("./data/");
	OutDirIO = std::string(argv[optind +1]);

	/*-----random number initialization-----*/
	uint seedIO = 0;
	seedIO = atoi(argv[optind + 2]);

  /*---- Flux data ----------------*/
  // Already loaded


	/*-----Geneartion-----*/

	VectGenIO *io = new VectGenIO(OutDirIO, seedIO, FluxFile, nuEne_min, nuEne_max);

	io->DoProcess(NumEv);

}
