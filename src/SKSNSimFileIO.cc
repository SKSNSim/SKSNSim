/**********************************
 * File: SKSNSimFileIO.cc
 **********************************/

#include <set>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <TFile.h>
#include <TFileCacheWrite.h>
#include <TTree.h>
#include <TMath.h>
#include <mcinfo.h>
#include "SKSNSimFileIO.hh"
#include "skrun.h"
#include "SKSNSimTools.hh"

/*
 * PROTOTYPE DECLARARATION for internal use
 */
std::vector<std::tuple<int,int,int>> ReadTimeEventFile(TRandom &rndgen, int runBegin, int runEnd, double nev_per_min = 1.0);

void SKSNSimFileOutTFile::Open(const std::string fname, const bool including_snevtinfo) {
  if(m_fileptr != NULL){
    Close();
    delete m_fileptr;
  }

  m_fileptr = new TFile(fname.c_str(), "RECREATE");
  std::cout << "Opened : " << fname << " (" << m_fileptr << ") "<< std::endl;
	//------------------------------------------------------------------------
	// set write cache to 40MB 
	TFileCacheWrite *cachew = new TFileCacheWrite(m_fileptr, 40*1024*1024);
	//------------------------------------------------------------------------

	/*-----define class-----*/
	m_MC = new MCInfo;
	m_MC->Clear();
	m_MC->SetName("MC");

  if(including_snevtinfo){
    m_SN = new SNEvtInfo;
    m_SN->Clear();
    m_SN->SetName("SN");
  }

	/*-----define branch-----*/
	TList *TopBranch = new TList;
	TopBranch->Add(m_MC);
  if(including_snevtinfo) TopBranch->Add(m_SN);

	// define tree
	m_OutTree = new TTree("data", "SK 5 tree");
	m_OutWeightTree = new TTree("weightTr", "Weight");
	// new MF
	m_OutTree->SetCacheSize(40*1024*1024);
	int bufsize = 8*1024*1024;      // may be this is the best 15-OCT-2007 Y.T.
	m_OutTree->Branch(TopBranch,bufsize);
  m_OutWeightTree->Branch("weight", &weight, "weight/D");
}

void SKSNSimFileOutTFile::Close(){
  std::cout << m_fileptr->GetName() << std::endl;
	m_fileptr->cd();
	m_OutTree->Write();
  m_OutWeightTree->Write();
  m_fileptr->Save();
	m_OutTree->Reset();
  m_OutWeightTree->Reset();
	m_fileptr->Close();
}

void SKSNSimFileOutTFile::Write(const SKSNSimSNEventVector &ev){

  if(ev.GetSNEvtInfoIEvt() != -1 && m_SN != NULL){
    m_SN->iEvt  = ev.GetSNEvtInfoIEvt();
    m_SN->rType = ev.GetSNEvtInfoRType();
    m_SN->nuType = ev.GetSNEvtInfoNuType();
    m_SN->rTime = ev.GetSNEvtInfoRTime();
    m_SN->nuEne = ev.GetSNEvtInfoNuEne();
    m_SN->nuDir[0] = ev.GetSNEvtInfoNuDir(0);
    m_SN->nuDir[1] = ev.GetSNEvtInfoNuDir(1);
    m_SN->nuDir[2] = ev.GetSNEvtInfoNuDir(2);
    m_SN->rVtx[0] = ev.GetSNEvtInfoRVtx(0);
    m_SN->rVtx[1] = ev.GetSNEvtInfoRVtx(1);
    m_SN->rVtx[2] = ev.GetSNEvtInfoRVtx(2);
  }

  m_MC->mcrun = ev.GetRunnum();
  m_MC->mcninfo = 2;
  m_MC->mcinfo[1] = ev.GetSubRunnum();

  // MCVERTEX (see $SKOFL_ROOT/inc/vcvrtx.h )
  m_MC->nvtxvc = ev.GetNVertex();
  for(int i = 0; i < ev.GetNVertex(); i++){
    m_MC->pvtxvc[i][0] = ev.GetVertexPositionX(i);
    m_MC->pvtxvc[i][1] = ev.GetVertexPositionY(i);
    m_MC->pvtxvc[i][2] = ev.GetVertexPositionZ(i);
    m_MC->timvvc[i] = ev.GetVertexTime(i);
    m_MC->iflvvc[i] = ev.GetVertexIFLVVC(i);
    m_MC->iparvc[i] = ev.GetVertexIPARVC(i);
  }


  m_MC->nvc = ev.GetNTrack();
  for(int i = 0; i < ev.GetNTrack(); i++){
    m_MC->ipvc[i] = ev.GetTrackPID(i);
    m_MC->energy[i] = ev.GetTrackEnergy(i);
    m_MC->pvc[i][0] = ev.GetTrackMomentumX(i);
    m_MC->pvc[i][1] = ev.GetTrackMomentumY(i);
    m_MC->pvc[i][2] = ev.GetTrackMomentumZ(i);
    m_MC->iorgvc[i] = ev.GetTrackIORGVC(i);
    m_MC->ivtivc[i] = ev.GetTrackIVTIVC(i);
    m_MC->ivtfvc[i] = ev.GetTrackIVTFVC(i);
    m_MC->iflgvc[i] = ev.GetTrackIFLGVC(i);
    m_MC->icrnvc[i] = ev.GetTrackICRNVC(i);
  }

  m_OutTree->Fill();

  weight = ev.GetWeight();
  m_OutWeightTree->Fill();
}

std::vector<SKSNSimFileSet> GenerateOutputFileListNoRuntime(const SKSNSimUserConfiguration &conf){
  std::vector<SKSNSimFileSet> flist;
  if( conf.GetNormRuntime() ) return flist;
  int nfile = conf.GetNumEvents() / conf.GetNumEventsPerFile();
  int nfile_rem = conf.GetNumEvents() % conf.GetNumEventsPerFile();

  const std::string suffix = ( conf.GetOFileMode() == SKSNSimUserConfiguration::MODEOFILE::kSKROOT ? "root" : "txt");

  auto genFileName = [](const SKSNSimUserConfiguration &c, int i, const std::string suffix){
    char buf[1000];
    snprintf(buf, 999, "%s/%s_%06d.%s", c.GetOutputDirectory().c_str(), c.GetOutputPrefix().c_str(), i, suffix.c_str());
    if( c.GetOutputNameTemplate() != SKSNSimUserConfiguration::GetDefaultOutputNameTemplate() ){
      auto pos = c.GetOutputNameTemplate().find("RUNNUM");
      auto pref = c.GetOutputNameTemplate().substr(0,pos);
      auto suf  = c.GetOutputNameTemplate().substr(pos+6);
      snprintf(buf, 999, "%s/%s%06d%s", c.GetOutputDirectory().c_str(), pref.c_str(), i, suf.c_str());
    }
    return std::string(buf);
  };

  for(int i = 0; i < nfile; i++){
    auto fn = genFileName(conf, i, suffix);
    flist.push_back(SKSNSimFileSet(fn, conf.GetRunnum(), conf.GetSubRunnum(), conf.GetNumEventsPerFile()));
  }
  if( nfile_rem != 0){
    auto fn = genFileName(conf, nfile, suffix);
    flist.push_back(SKSNSimFileSet(fn, conf.GetRunnum(), conf.GetSubRunnum(), nfile_rem));
  }
  return flist;
}

std::vector<SKSNSimFileSet> GenerateOutputFileListRuntime(SKSNSimUserConfiguration &conf){
  std::vector<SKSNSimFileSet> buffer;
  if( !conf.GetNormRuntime() ) return buffer;

  std::vector<std::tuple<int,double>> livetimes;
  if( conf.CheckMODERuntime() == SKSNSimUserConfiguration::MODERUNTIME::kRUNTIMERUNNUM ){
    livetimes = SKSNSimLiveTime::LoadLiveTime( conf.GetRuntimeRunBegin(), conf.GetRuntimeRunEnd());
  } else if( conf.CheckMODERuntime() == SKSNSimUserConfiguration::MODERUNTIME::kRUNTIMEPERIOD ){
    livetimes = SKSNSimLiveTime::LoadLiveTime( (SKSNSIMENUM::SKPERIOD) conf.GetRuntimePeriod());
  }
  std::vector<std::tuple<int,int>> nev_livetimes = SKSNSimLiveTime::ConvertExpectedEvt(*conf.GetRandomGenerator(), livetimes, conf.GetRuntimeFactor());

  for(auto it = nev_livetimes.begin(); it != nev_livetimes.end(); it++){
    auto vectio = std::make_unique<SKSNSimFileOutTFile>();
    const int run = std::get<0>(*it);
    std::stringstream ss;
    if( conf.GetOutputNameTemplate() != SKSNSimUserConfiguration::GetDefaultOutputNameTemplate() ){
      auto pos  = conf.GetOutputNameTemplate().find("RUNNUM");
      auto pref = conf.GetOutputNameTemplate().substr(0,pos);
      auto suf  = conf.GetOutputNameTemplate().substr(pos+6);
      ss << conf.GetOutputDirectory() << "/" << pref << std::setfill('0') << std::setw(6) << run << suf;

    } else {
      std::string suffix = (conf.GetOFileMode() == SKSNSimUserConfiguration::MODEOFILE::kSKROOT ? "root" : "txt");
      ss << conf.GetOutputDirectory() << "/" << conf.GetOutputPrefix() << "." << std::setfill('0') << std::setw(6) << run << "." << suffix;
    }
    const std::string fname(ss.str());

    std::cout << "fname " << fname << " nev " << std::get<1>(*it) <<std::endl;

    int subrun = 0;
    int nev = std::get<1>(*it);
    buffer.push_back(SKSNSimFileSet(fname, run,subrun, nev));
  }


  /* Old TimeEvent loader 
  const int runBegin = conf.GetRuntimeRunBegin();
  const int runEnd   = conf.GetRuntimeRunEnd();

  auto runev = ReadTimeEventFile(*conf.GetRandomGenerator(), runBegin, runEnd, conf.GetRuntimeFactor());

  auto genFileName = [](const SKSNSimUserConfiguration &c, int r,int sr){
    char buf[1000];
    snprintf(buf, 999, "%s/%s_r%06d_%06d.root", c.GetOutputDirectory().c_str(), c.GetOutputPrefix().c_str(), r, sr);
    return std::string(buf);
  };
  for(auto it = runev.begin(); it != runev.end(); it++){
    int run = std::get<0>(*it);
    int subrun = std::get<1>(*it);
    int nev = std::get<2>(*it);
    auto fn = genFileName(conf, run, subrun);
    buffer.push_back(SKSNSimFileSet(fn, run,subrun, nev));
  }
  */
  return buffer;
}

std::vector<SKSNSimFileSet> GenerateOutputFileList(SKSNSimUserConfiguration &conf){
  std::vector<SKSNSimFileSet> flist;
  if( conf.GetNormRuntime() )
    flist = GenerateOutputFileListRuntime(conf);
  else
    flist = GenerateOutputFileListNoRuntime(conf);

  return flist;
}

SKRUN FindSKPeriod(int run){
    if ( run < SK_IV_BEGIN )
      return SK_I_II_III_BEGIN;
    else if ( SK_IV_BEGIN <= run && run <= SK_IV_END )
      return SK_IV_BEGIN;
    else if ( SK_V_BEGIN <= run && run <= SK_V_END ) 
      return SK_V_BEGIN;
    else if ( SK_VI_BEGIN <= run && run <= SK_VI_END)
      return SK_VI_BEGIN;
    else if (SK_VII_BEGIN <= run )
      return SK_VII_BEGIN;
    
    return SK_I_II_III_BEGIN;
}

std::string FindTimeFile(int run) {
    if ( run < SK_IV_BEGIN || run >= SK_VII_BEGIN ) {
      std::cout << "reference run number is not correct"<<std::endl;
      return "";
    } 
    else if ( SK_IV_BEGIN <= run && run < SK_IV_END ) { 
      return "/home/sklowe/realtime_sk4_rep/solar_apr19/timevent/livesubruns.r061525.r077958"; 
    }
    else if ( SK_V_BEGIN <= run && run < SK_V_END ) { 
      return "/home/sklowe/realtime_sk5_rep/solar_nov20/timevent/livesubruns.r080539.r082915";
    }
    else if ( SK_VI_BEGIN <= run && run < SK_VI_END ) { 
      return "/home/sklowe/realtime_sk6_rep/solar_may22/timevent/livesubruns.r085000.r087073"; 
    }
    return "";
}

void LoadTimeEventsFromFile(std::vector<std::tuple<int,int,double>> &buffer, SKRUN p){
  const std::string fname = FindTimeFile(p);
  buffer.clear();
  if(fname=="") return;

  std::ifstream ifs(fname);
  if(!ifs.is_open()){
    std::cerr << "Failed to open : " << fname << std::endl;
    return;
  }

  std::cout << "File opend: " << fname << std::endl;
  for(std::array<char,500> b; ifs.getline(&b[0], 500); ){
    int r, sr;
    double t, tt;
    sscanf(&b[0], "%d %d %lf %lf", &r, &sr, &t, &tt);
    buffer.push_back( std::make_tuple(r,sr,t));
  }

  ifs.close();
  return;
}

std::vector<std::tuple<int,int,int>> ReadTimeEventFile(TRandom &rndgen, int runBegin, int runEnd, double nev_per_min)
{
  std::cout <<" Estimate # of event from timevent file "<<std::endl;
  std::vector<std::tuple<int,int,int>> runev; // run,subrun,nev
  std::vector<std::tuple<int,int,double>> runtime; // run,subrun,time

  auto skperiodBegin = FindSKPeriod(runBegin);
  auto skperiodEnd   = FindSKPeriod(runEnd);
  std::cout << "Period run = [" << skperiodBegin << ", " << skperiodEnd << ")" << std::endl;
  std::set<SKRUN> skperiods{SK_IV_BEGIN, SK_V_BEGIN, SK_VI_BEGIN, SK_VII_BEGIN};
  for(auto it = skperiods.begin(); it != skperiods.end();){
    if(skperiodBegin <= *it && *it <= skperiodEnd ) it++;
    else 
      it = skperiods.erase(it);
  }
  
  for(auto it = skperiods.begin(); it != skperiods.end(); it++)
    LoadTimeEventsFromFile( runtime, *it);

  for(auto it = runtime.begin(); it != runtime.end(); ){
    if( runBegin <= std::get<0>(*it) && std::get<0>(*it) < runEnd ) it++;
    else it = runtime.erase(it);
  }

  auto convDoubleToInt = [] (TRandom &r, double d, double w){
    double nev_double = d * w;
    int nev = int(nev_double);
    double nev_frac = nev_double - nev; 
    double rand = r.Rndm();
    int nev_frac_int = 0;
    if(nev_frac < rand) nev_frac_int = 0;
    else nev_frac_int = 1;
    return nev+nev_frac_int;
  };

  const double weight = nev_per_min / 60.0;
  for(auto it = runtime.begin(); it != runtime.end(); it++)
    runev.push_back(
        std::make_tuple( std::get<0>(*it), std::get<1>(*it),
          convDoubleToInt(rndgen, std::get<2>(*it), weight)));

  return runev;
}

void SKSNSimFileOutNuance::Open ( const std::string fname) {
  ofs.reset(new std::ofstream( fname ));
  if( !ofs->good() ) { 
    std::cout << "ERRROR: failed to open output file ( " << fname << std::endl;
    return;
  }
  std::cout << "Opened output file: " << fname << std::endl;
}

void SKSNSimFileOutNuance::Close() {
  *ofs << "stop" << std::endl;
  ofs->close();
}

void SKSNSimFileOutNuance::Write(const SKSNSimSNEventVector &ev) {

  auto convReactionMode = [](const SKSNSimSNEventVector &ev) {
    return "nuance " + std::to_string(ev.GetSNEvtInfoRType());};
  auto convVertex = [] (const SKSNSimSNEventVector &ev) {
    return ( "vertex "
        //+ std::to_string(ev.GetSNEvtInfoRVtx(0)) + " "
        //+ std::to_string(ev.GetSNEvtInfoRVtx(1)) + " "
        //+ std::to_string(ev.GetSNEvtInfoRVtx(2)) + " " 
        + std::to_string(ev.GetVertexPositionX(0)) + " "
        + std::to_string(ev.GetVertexPositionY(0)) + " "
        + std::to_string(ev.GetVertexPositionZ(0)) + " "
        + std::to_string(ev.GetSNEvtInfoRTime())
          );
  };
  auto convTrack = [] (const SKSNSimSNEventVector &ev, int i) {
    double dir[3];
    double totmon = ev.GetTrackMomentumX(i)*ev.GetTrackMomentumX(i);
    totmon += ev.GetTrackMomentumY(i)*ev.GetTrackMomentumY(i);
    totmon += ev.GetTrackMomentumZ(i)*ev.GetTrackMomentumZ(i);
    totmon = TMath::Sqrt(totmon);
    if( totmon <= 0.0){
      dir[0] = 0.0; dir[1] = 0.0; dir[2] = 0.0;
    } else {
      dir[0] = ev.GetTrackMomentumX(i) / totmon;
      dir[1] = ev.GetTrackMomentumY(i) / totmon;
      dir[2] = ev.GetTrackMomentumZ(i) / totmon;
    }
    int mode = -9999;
    if( ev.GetTrackICRNVC(i) != 0 ) mode = 0;
    else if ( ev.GetTrackIORGVC(i) != 0 ) mode = -1;
    else mode = -2;
    return "track "
      + std::to_string(ev.GetTrackPID(i)) + " "
      + std::to_string(ev.GetTrackEnergy(i)) + " "
      + std::to_string(dir[0]) + " " + std::to_string(dir[1]) + " " + std::to_string(dir[2]) + " "
      + std::to_string(mode);
  };


  *ofs << "begin" << std::endl
    << convReactionMode(ev) << std::endl
    << convVertex(ev) << std::endl;
  for(int i = 0; i < ev.GetNTrack(); i++) *ofs << convTrack(ev, i) << std::endl;
  *ofs << "end" << std::endl;

}

