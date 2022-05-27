
#include "VectGenUtil.hh"
#include "VectGenSetBin.hh"
#include "VectGenIO.hh"

VectGenIO::VectGenIO(std::string ModelName, int nuosc, double distIO, int flagIO, std::string OutDirIO, uint seedIO)
{

	RatioTo10kpc = 1. / SQ(distIO);
	flag_event = flagIO;
	generator = new TRandom3(seedIO);
	/*
	   if( seedIO == 0 ){ // generating the seed
	   seedIO = generator->GetSeed();
	//std::cout << seedIO << std::endl;
	generator->SetSeed(seedIO);
	}
	*/

	/*-----determine SN direction-----*/
	int date[3], time[3];
	date[0] = 2011;
	date[1] = 3;
	date[2] = 23;
	time[0] = 0;
	time[1] = 0;
	time[2] = 0;
	float sdir[3], ra, dec;
	sn_sundir_( date, time, sdir, & ra, & dec);
	for( int i = 0; i < 3; i++ ){ snDir[i] = ( double )sdir[i]; }

	double theta = acos( snDir[2] );
	double phi = atan2( snDir[1],  snDir[0] );

	Rmat[0][0] = cos( theta ) * cos( phi );
	Rmat[0][1] = -sin( phi );
	Rmat[0][2] = sin( theta ) * cos( phi );
	Rmat[1][0] = cos( theta ) * sin( phi );
	Rmat[1][1] = cos( phi );
	Rmat[1][2] = sin( theta ) * sin( phi );
	Rmat[2][0] = -sin( theta );
	Rmat[2][1] = 0.;
	Rmat[2][2] = cos( theta );

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
	ssname1 << "/home/sklowe/supernova/data/" << ModelName;
	std::string FileIn = ssname1.str();
	nuflux = new VectGenSnNakazato(FileIn.c_str());
	//nuflux = new SnWilson();

	//Class call of cross-section calculation
	nucrs = new VectGenNuCrosssection();
	nucrs->ReadCsNuElastic();

	//Class call of cross-section calculation
	ocrs = new VectGenOxyCrosssection();
	osub = new VectGenOxyCrosssectionSub();


	// set text file name
	std::stringstream ssname3;
	//ssname3 << OutDirIO << "/" << ModelName << ".dat";
	//ssname3 << OutDirIO << "/" << ModelName << "/eventrate.dat";
	ssname3 << OutDirIO << "/eventrate.dat";
	FileText = ssname3.str();
	std::cout << "text file name " << FileText << std::endl;

	if(flag_event == 1){ // for making events

	  // set output directory for vector root file
	  char del ='/';
	  auto subStr = split(ModelName, del);
	  std::stringstream ssname4;
	  //ssname4 << OutDirIO << subStr[0] << "/event/";
	  //ssname4 << OutDirIO << "/" << ModelName << "/event/";
	  ssname4 << OutDirIO << "/event/";
	  OutDir = ssname4.str();
	  std::cout << OutDir << std::endl;
	  
	}

}

VectGenIO::VectGenIO(std::string OutDirIO, uint seedIO, std::string FluxFileName)
{
	generator = new TRandom3(seedIO);
	/*
	   if( seedIO == 0 ){ // generating the seed
	   seedIO = generator->GetSeed();
	//std::cout << seedIO << std::endl;
	generator->SetSeed(seedIO);
	}
	*/

	/*-----set binning-----*/
	VectGenSetBinValues();

	//Class call of flux/spectrum table
	//(should be here, but now in VectGenGenerator::Process(int)
  nuflux_dsbn.reset(new FluxCalculation(FluxFileName));

	//Class call of cross-section calculation
	nucrs = new VectGenNuCrosssection();
	nucrs->ReadCsNuElastic();

	// set output directory for vector root file
	/*
	char del ='/';
	auto subStr = split(ModelName, del);
	std::stringstream ssname4;
	//ssname4 << OutDirIO << subStr[0] << "/event/";
	//ssname4 << OutDirIO << "/" << ModelName << "/event/";
	ssname4 << OutDirIO << "/event/";
	OutDir = ssname4.str();
	*/
	OutDir = OutDirIO;
	std::cout << OutDir << std::endl;
}

void VectGenIO::DoProcess()
{

	Process();

	//std::cout << "each event number is "<< totNuebarp << " " << totNueElastic << " " << totNuebarElastic << " " << totNuxElastic << " " << totNuxbarElastic << std::endl;

	std::ofstream ofs(FileText.c_str());
	ofs << totNuebarp << " " << totNueElastic << " " << totNuebarElastic << " " << totNuxElastic << " " << totNuxbarElastic << std::endl;

}

void VectGenIO::DoProcess(int NumEv)
{
	Process(NumEv);
}
