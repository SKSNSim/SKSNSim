
#include "VectGenIO.hh"
#include <getopt.h>
#include <string>
#include <iostream>

constexpr double DEFAULT_NUENE_MIN = 5.;
constexpr double DEFAULT_NUENE_MAX = 100.;

void ShowUsage(char *arg0){
  std::cerr<< arg0 << " [options] {num_eve} {out_dir} {seed}" << std::endl;
  std::cerr<< "Usage: vector generetaor for DSNB MC" << std::endl;
  std::cerr<<" Argument for "<< arg0 << std::endl;
  //std::cerr<<" 1st: Number of event" << std::endl;
  //std::cerr<<" 2nd: Output directory, default: ./data/"<< std::endl;
  //std::cerr<<" 3rd: Random seed"<< std::endl;
  std::cerr<<" Options: -f(--flux_file), -h(--help), -r(--ref_runnum), --nu_ene_min, --nu_ene_max \n"
    << "\t-f(--flux_file) [flux_data_file]: specify flux data file to be used" << std::endl
    << "\t-r(--ref_runnum) [reference_run_number]: specify run number stored to output file" << std::endl
    << "\t-s(--seed) [random seed]: (mandatory)specify random seed number" << std::endl
    << "\t-n(--num_event) [number of event]: specify the number of events generated" << std::endl
    << "\t--out [output_vector_name]: (mandatory) specify the name for output vector file" << std::endl
    << "\t--nu_ene_min [energy_in_MeV]: minimum total energy [MeV] to be generated (default = " << DEFAULT_NUENE_MIN << " MeV)" << std::endl
    << "\t--TIMEEVENT [true/false]: if the number of event is decided from live time" << std::endl;
  std::cerr<<" return -1"<<std::endl;
}

int main( int argc, char ** argv )
{

	std::cout<<" Start Generating Event "<< argc <<std::endl;

  //=====================
  // argument controll
  int c, digit_optind = 0;
  double nuEne_min = -9999.;
  double nuEne_max = 9999.;
  std::string FluxFile ("dsnb_flux/horiuchi/8MeV_Nominal.dat");
	uint seed = 0;
  int numEvent = 0;
  int RefRunNum = 0;
  bool useTimeEvent = false;
	std::string OutputFile("");
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    unsigned long seed_read = 0;
    enum OPTIONS {kNuEneMin = 0, kNuEneMax, kFluxFile, kRefRunNum, kSeed, kNumEvent, kOutFile, kTimeEvent, kHelp};
    static struct option long_options[] = {
      {"nu_ene_min",  required_argument, 0, 0},
      {"nu_ene_max",  required_argument, 0, 0},
      {"flux_file",   required_argument, 0, 0},
      {"ref_runnum",  required_argument, 0, 0},
      {"seed",        required_argument, 0, 0},
      {"num_event",   required_argument, 0, 0},
      {"out",         required_argument, 0, 0},
      {"TIMEEVENT",   required_argument, 0, 0},
      {"help",        no_argument, 0, 0},
      {0,             0, 0, 0}
    };

    c = getopt_long( argc, argv, "hf:r:s:n:v",
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

          case kRefRunNum:
            RefRunNum = std::stoi(optarg);
            std::cerr << "GETOPT: reference run number = " << RefRunNum << std::endl;
            break;
            
          case kSeed:
            seed = std::stoi(optarg);
            std::cerr << "GETOPT: random seed = " << seed << std::endl;
            break;
            
          case kNumEvent:
            numEvent = std::stoi(optarg);
            std::cerr << "GETOPT: number of generation event = " << numEvent << std::endl;
            break;
            
          case kOutFile:
            OutputFile = std::string(optarg);
            std::cerr << "GETOPT: output file name = " << OutputFile << std::endl;
            break;

          case kTimeEvent:
            useTimeEvent = std::stoi(optarg);
            std::cerr << "GETOPT: Use livetime info. for number of event = " << useTimeEvent << std::endl;
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

      case 'r':
        RefRunNum = std::stoi(optarg);
        std::cerr << "GETOPT: reference run number = " << RefRunNum << std::endl;
        break;

      case 's':
        seed = std::stoi(optarg);
        std::cerr << "GETOPT: random seed = " << seed << std::endl;
        break;

      case 'n':
        numEvent = std::stoi(optarg);
        std::cerr << "GETOPT: number of generation event = " << numEvent << std::endl;
        break;
            
      case 'h':
        ShowUsage(argv[0]);
        return EXIT_SUCCESS;
        break;

      default:
        std::cerr << "GETOPT: strange argument: code = " << c << " is ignored !!"<<std::endl;
    }
  }

//  if( optind+3 > argc){
//    ShowUsage(argv[0]);
//    return EXIT_FAILURE;
//  }

  if ( OutputFile == "") {
    std::cerr <<" Output file name must be specified"<<std::endl; 
    ShowUsage(argv[0]);
    return EXIT_FAILURE;
  }



//	/*-----Number of Generated Event-----*/
//	int NumEv = atoi(argv[optind]);

//	/*-----Output directory-----*/
//	std::string OutDirIO("./data/");
//	OutDirIO = std::string(argv[optind +1]);

//	/*-----random number initialization-----*/
//	uint seedIO = 0;
//	seedIO = atoi(argv[optind + 2]);

  /*---- Flux data ----------------*/
  // Already loaded


	/*-----Geneartion-----*/

//	VectGenIO *io = new VectGenIO(OutDirIO, seedIO, FluxFile);
	VectGenIO *io = new VectGenIO(seed);
  io->SetFluxFile(FluxFile);
  io->OpenOutputFile(OutputFile);
  io->SetRefRunNumber(RefRunNum);
  io->SetUseTimeEvent(useTimeEvent);

	io->DoProcess(numEvent);

  io->CloseOutputFile();

}
