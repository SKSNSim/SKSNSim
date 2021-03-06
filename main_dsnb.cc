
#include "VectGenIO.hh"
#include "VectGenSetBin.hh"
#include <getopt.h>
#include <string>
#include <iostream>
#include <sstream>

const double DEFAULT_NUENE_MIN = nuEneMin;
const double DEFAULT_NUENE_MAX = nuEneMax;
const std::string DEFAULT_FLUX_FILE ("dsnb_flux/horiuchi/8MeV_Nominal.dat");

void ShowUsage(char *arg0){
  std::cerr<< arg0 << " [options] --out {out_dir} -s {seed}" << std::endl;
  std::cerr<< "Usage: vector generetaor for DSNB MC" << std::endl;
  std::cerr<<" Argument for "<< arg0 << std::endl;
  //std::cerr<<" 1st: Number of event" << std::endl;
  //std::cerr<<" 2nd: Output directory, default: ./data/"<< std::endl;
  //std::cerr<<" 3rd: Random seed"<< std::endl;
  std::cerr<<" Options: -f(--flux_file), -h(--help), -r(--ref_runnum), -o(--out), -n, --NEVPERFILE, --nu_ene_min, --nu_ene_max \n"
    << "\t-f(--flux_file) [flux_data_file]: specify flux data file to be used" << std::endl
    << "\t-r(--ref_runnum) [reference_run_number]: specify run number stored to output file" << std::endl
    << "\t-s(--seed) [random seed]: (mandatory)specify random seed number" << std::endl
    << "\t-n(--num_event) [number of event]: specify the number of events generated" << std::endl
    << "\t--NEVPERFILE [number of event]: number of events per one file (default: all events are stored in single file" << std::endl
    << "\t-o(--out) [output_vector_name]: (mandatory) specify the name for output vector file" << std::endl
    << "\t--nu_ene_min [energy_in_MeV]: minimum total energy [MeV] to be generated (default = " << DEFAULT_NUENE_MIN << " MeV)" << std::endl
    << "\t--nu_ene_max [energy_in_MeV]: maximum total energy [MeV] to be generated (default = " << DEFAULT_NUENE_MAX << " MeV)" << std::endl
    << "\t--TIMEEVENT [true(1)/false(0)]: if the number of event is decided from live time" << std::endl
    << "\t--FLATFLUX  [true(1)/false(0)]: if the flux assume flat " << std::endl;
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
  std::string FluxFile (DEFAULT_FLUX_FILE);
	uint seed = 0;
  int numEvent = 0;
  int numEventPerFile = 0;
  int RefRunNum = 0;
  bool useTimeEvent = false;
  bool useFlatFlux = false;
	std::string OutputFile("");
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    enum OPTIONS {kNuEneMin = 0, kNuEneMax, kFluxFile, kRefRunNum, kSeed, kNumEvent, kNumEventPerFile, kOutFile, kTimeEvent, kFlatFlux, kHelp}; // this order should be same with long_options below
    static struct option long_options[] = {
      {"nu_ene_min",  required_argument, 0, 0},
      {"nu_ene_max",  required_argument, 0, 0},
      {"flux_file",   required_argument, 0, 0},
      {"ref_runnum",  required_argument, 0, 0},
      {"seed",        required_argument, 0, 0},
      {"num_event",   required_argument, 0, 0},
      {"NEVPERFILE",  required_argument, 0, 0},
      {"out",         required_argument, 0, 0},
      {"TIMEEVENT",   required_argument, 0, 0},
      {"FLATFLUX",    required_argument, 0, 0},
      {"help",        no_argument, 0, 0},
      {0,             0, 0, 0}
    };

    c = getopt_long( argc, argv, "hf:r:s:n:o:v",
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

          case kNumEventPerFile:
            numEventPerFile = std::stoi(optarg);
            std::cerr << "GETOPT: number of event per file = " << numEventPerFile << std::endl;
            break;

          case kOutFile:
            OutputFile = std::string(optarg);
            std::cerr << "GETOPT: output file name = " << OutputFile << std::endl;
            break;

          case kTimeEvent:
            useTimeEvent = std::stoi(optarg);
            std::cerr << "GETOPT: Use livetime info. for number of event = " << useTimeEvent << std::endl;
            break;

          case kFlatFlux:
            useFlatFlux= std::stoi(optarg);
            std::cerr << "GETOPT: Use flat flux for positron spectrum = " << useFlatFlux << std::endl;
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

      case 'o': // same with kOutFile
        OutputFile = std::string(optarg);
        std::cerr << "GETOPT: output file name = " << OutputFile << std::endl;
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


	/*-----Geneartion-----*/
	VectGenIO *io = new VectGenIO(seed);
  if (nuEne_min>-9998. && nuEne_max< 9998.) io->SetNuEnergyRange(nuEne_min, nuEne_max);
  else if (nuEne_min>-9998. && nuEne_max> 9998.)
  {
    std::cerr << "[error] spcify upper limit also by --nu_ene_max {x_MeV}" << std::endl;
    return EXIT_FAILURE;
  }
  else if (nuEne_min<-9998. && nuEne_max< 9998.)
  {
    std::cerr << "[error] spcify lower limit also by --nu_ene_min {x_MeV}" << std::endl;
    return EXIT_FAILURE;
  }
  io->SetFluxFile(FluxFile);

  io->SetRefRunNumber(RefRunNum);
  io->SetUseTimeEvent(useTimeEvent);
  io->SetUseFlatFlux(useFlatFlux);

  if( numEventPerFile == 0 ) { // push all events in one file: OutputFile
    io->OpenOutputFile(OutputFile);
    io->DoProcess(numEvent);
    io->CloseOutputFile();
  } else {
    const int nfile = numEvent / numEventPerFile;
    const int nEventRem = numEvent % numEventPerFile;
    for( int f = 0; f < nfile; f++) {
      std::stringstream ss;
      ss << OutputFile << "/" << std::setw(6) << std::setfill('0') << f << ".root";
      std::cout << "Output to " << ss.str() << std::endl;
      io->OpenOutputFile(ss.str());
      io->DoProcess(numEventPerFile);
      io->CloseOutputFile();
    }
    if ( nEventRem != 0) {
      std::stringstream ss;
      ss << OutputFile << "/" << std::setw(6) << std::setfill('0') << nfile << ".root";
      std::cout << "Output to " << ss.str() << std::endl;
      io->OpenOutputFile(ss.str());
      io->DoProcess(nEventRem);
      io->CloseOutputFile();
    }
  }

}
