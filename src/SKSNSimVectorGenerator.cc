/***********************************
 * SKSNSimVectorGenerator.cc
 * *********************************/
#include <functional>
#include <algorithm>
#include <TRandom3.h>
#include <geotnkC.h>
#include "SKSNSimVectorGenerator.hh"
#include "SKSNSimConstant.hh"
#include "SKSNSimCrosssection.hh"
#include "SKSNSimTools.hh"

using namespace SKSNSimPhysConst;


double SKSNSimVectorGenerator::FindMaxProb ( SKSNSimFluxModel &flux, SKSNSimCrosssectionModel &xsec){
  double maxP = 0.;
  const double ene_min = flux.GetEnergyLimitMin();
  const double ene_max = flux.GetEnergyLimitMax();
  constexpr double cost_min = -1.;
  constexpr double cost_max =  1.;
  constexpr size_t nbin_ene = 1000; // TODO define appropriate binsize
  constexpr size_t nbin_cost = 1000; // TODO define appropriate binsize
  const double diff_ene = (ene_max - ene_min)/nbin_ene;
  constexpr double diff_cost = (cost_max - cost_min)/nbin_cost;
  for(size_t i = 0; i < nbin_ene; i++){
    const double ene =  ene_min +  diff_ene* double(i);
    for(size_t j = 0; j < nbin_cost; j++){
      const double cost =  cost_min +  diff_cost* double(i);
      const double p = flux.GetFlux(ene, 0.0, SKSNSimFluxModel::FLUXNUEB) * xsec.GetDiffCrosssection(ene, cost).first;
      if(maxP < p) maxP = p;
    }
  }
  return maxP;
}


double SKSNSimVectorGenerator::SetMaximumHitProbability(){
  // TODO fix to use appropriate combination of flux and xsec models
  if( fluxmodels.size() < 1 || xsecmodels.size() < 1 ){
    m_max_hit_probability = -1.0;
    return m_max_hit_probability;
  }
  m_max_hit_probability = FindMaxProb(*fluxmodels[0], *xsecmodels[0]);
  return m_max_hit_probability;
}

SKSNSimSNEventVector SKSNSimVectorGenerator::GenerateEventIBD() {
  SKSNSimSNEventVector ev;

  double nuEne, cost, eEne;
  double nEne;

  TRandom &rng = *randomgenerator;
  if(fluxmodels.size() == 0) return ev;
  SKSNSimFluxModel &flux = *fluxmodels[0]; // TODO modify for user to select models
  if(xsecmodels.size() == 0) return ev;
  SKSNSimCrosssectionModel &xsec = *xsecmodels[0]; // TODO modify for user to select models

  auto SQ = [](double a){ return a*a;};

  // determine neutrino and positron energy, and its direction
  while( 1 ){
    nuEne = rng.Uniform( GetEnergyMin(), GetEnergyMax());

    const double nuFlux = flux.GetFlux(nuEne, 0.0, SKSNSimFluxModel::FLUXNUEB);

    double sigm;
    cost = rng.Uniform( -1., 1.);
    auto xsecpair = xsec.GetDiffCrosssection(nuEne, cost);
    eEne = xsecpair.second;
    sigm = xsecpair.first;

    double p = nuFlux * sigm;
    double x = rng.Uniform( 0., m_max_hit_probability);
    if( x < p ) break;
  }
#ifdef DEBUG
  std::cout << "In GenerateEventIBD: eEne = " << eEne << std::endl;
#endif

  // determine neutrino direction
  //double nuDir[3];
  double theta = rng.Uniform(M_PI);
  double phi = rng.Uniform(2. * M_PI);
  UtilVector3<double> nuDir ( std::sin(theta)*std::cos(phi),
                                        sin(theta)*sin(phi),
                                                 cos(theta));


  // Rotation matrix of neutrino direction
  UtilMatrix3<double> Rmat( cos(theta)*cos(phi), -sin(phi), sin(theta)*cos(phi),
                            cos(theta)*sin(phi),  cos(phi), sin(theta)*sin(phi),
                                    -sin(theta),        0.,          cos(theta));

  // interaction point
  auto determinePosition = std::bind([](TRandom &rng, SKSNSIMENUM::TANKVOLUME t)
  {
    double rPositionRange = RINTK;
    double hPositionRange = ZPINTK;
    const double FVCUT = 200.;
    switch (t)
    {
      case SKSNSIMENUM::TANKVOLUME::kIDFV: //Fiducial volume
        rPositionRange = RINTK - FVCUT;
        hPositionRange = ZPINTK - FVCUT;
        break;
      case SKSNSIMENUM::TANKVOLUME::kIDFULL: //entire ID volume
        rPositionRange = RINTK;
        hPositionRange = ZPINTK;
        break;
      case SKSNSIMENUM::TANKVOLUME::kTANKFULL: //entire detector volume (including OD)
        rPositionRange = DITKTK;
        hPositionRange = ZPTKTK;
        break;
      default: //entire ID volume
        rPositionRange = RINTK;
        hPositionRange = ZPINTK;
    }
    //random inside the full tank (32.5kton)
    const double r2 = rng.Uniform(1.) * rPositionRange*rPositionRange ;
    const double r = std::sqrt( r2 );
    const double phi = rng.Uniform( 2. * M_PI);
    const double x = r * std::cos( phi );
    const double y = r * std::sin( phi );
    const double z = -hPositionRange + rng.Uniform( 2.*hPositionRange);

    return UtilVector3<double>(x,y,z);
  }, std::placeholders::_1, m_generator_volume);
  const auto xyz = determinePosition(*randomgenerator);

  // Fill into class
  // MCVERTEX (see $SKOFL_ROOT/inc/vcvrtx.h )                                                                               
  ev.AddVertex(xyz.x, xyz.y, xyz.z,
      1, 0, 0.);


  // IBD interaction

  // Original neutrino
  const auto nuMomentum = nuEne * nuDir.Unit();
  const auto nuebarid = ev.AddTrack(
      -12, nuEne,
      nuMomentum.x, nuMomentum.y, nuMomentum.z,
      0, 0, 1, -1, 0
      );

  // Original proton
  const auto protonid = ev.AddTrack(
      2212, Mp,
      0., 0., 0.,
      0, 0, 1, -1, 0
      );

  // Positron
  double amom = sqrt(SQ( eEne ) - SQ( Me ));
  double eTheta = acos( cost );
  double ePhi = rng.Uniform( -M_PI,  M_PI);

  // conversion the positron direction along the neutrino direction

  UtilVector3<double> origVec(
      sin( eTheta ) * cos( ePhi ),
      sin( eTheta ) * sin( ePhi ),
      cos( eTheta ));

  const UtilVector3<double> eDir = Rmat * origVec;

  //------------------------------------------------------------------------

  // Positron
  const auto positronMomentum = amom * eDir.Unit();
  const auto positronid = ev.AddTrack(
      -11, eEne,
      positronMomentum.x, positronMomentum.y, positronMomentum.z,
      1, 1, 1, 0, 1
      );

  // Neutron
  const UtilVector3<double> neutronMomentum = nuMomentum - positronMomentum;
  const auto neutronid = ev.AddTrack(
      2112,
      std::sqrt(neutronMomentum.Mag2() + Mn*Mn),
      neutronMomentum.x, neutronMomentum.y, neutronMomentum.z,
      1, 1, 1, 0, 1
      );

  ev.SetRunnum(999999); // TODO modify
  ev.SetSubRunnum(0);  // TODO

  return ev;
}

SKSNSimVectorSNGenerator::SKSNSimVectorSNGenerator():
  m_runnum((int)SKSNSIMENUM::SKPERIODRUN::SKMC ),
  m_subrunnum(0),
  m_generator_energy_min(0.0),
  m_generator_energy_max(300.0),
  m_nu_energy_nbins(3000),
  m_generator_time_min (0.0),
  m_generator_time_max (20.0),
  m_time_nbins(20000),
  m_fill_event(true),
  m_generator_volume( SKSNSIMENUM::TANKVOLUME::kIDFULL ),
  m_nuosc_type( SKSNSIMENUM::NEUTRINOOSCILLATION::kNONE ),
  m_distance_kpc(10.)
{
  m_sn_date[0] = 2011;
  m_sn_date[1] = 3;
  m_sn_date[2] = 23;
  m_sn_time[0] = 0;
  m_sn_time[1] = 0;
  m_sn_time[2] = 0;

  xsecmodels[XSECTYPE::mXSECIBD]       = std::make_unique<SKSNSimXSecIBDSV>();
  xsecmodels[XSECTYPE::mXSECELASTIC]   = std::make_unique<SKSNSimXSecNuElastic>();
  xsecmodels[XSECTYPE::mXSECOXYGEN]    = std::make_unique<SKSNSimXSecNuOxygen>();
  xsecmodels[XSECTYPE::mXSECOXYGENSUB] = std::make_unique<SKSNSimXSecNuOxygenSub>();
  xsecmodels[XSECTYPE::mXSECOXYGENNC]  = std::make_unique<SKSNSimXSecNuOxygenNC>();
}

std::vector<SKSNSimSNEventVector> SKSNSimVectorSNGenerator::GenerateEvents(){
  std::vector<SKSNSimSNEventVector> evt_buffer;
  SKSNSimBinnedFluxModel &flux = dynamic_cast<SKSNSimBinnedFluxModel&>(*fluxmodels[0]); // TODO selectable flux
  if(&flux == NULL) {
    std::cerr << "In GenerateEvents() no appropriate flux model (binned flux)" << std::endl;
    evt_buffer.clear();
    return evt_buffer;
  }

  SKSNSimXSecIBDSV       &xsecibd         = dynamic_cast<SKSNSimXSecIBDSV&>(      *xsecmodels[XSECTYPE::mXSECIBD]);
  SKSNSimXSecNuElastic   &xsecnuela       = dynamic_cast<SKSNSimXSecNuElastic&>(  *xsecmodels[XSECTYPE::mXSECELASTIC]);
  SKSNSimXSecNuOxygen    &xsecnuoxygen    = dynamic_cast<SKSNSimXSecNuOxygen&>(   *xsecmodels[XSECTYPE::mXSECOXYGEN]);
  SKSNSimXSecNuOxygenSub &xsecnuoxygensub = dynamic_cast<SKSNSimXSecNuOxygenSub&>(*xsecmodels[XSECTYPE::mXSECOXYGENSUB]);
  SKSNSimXSecNuOxygenNC  &xsecnuoxygennc  = dynamic_cast<SKSNSimXSecNuOxygenNC&>( *xsecmodels[XSECTYPE::mXSECOXYGENNC]);

	std::cout << "Prcess of sn_burst side" << std::endl;//nakanisi
	/*---- Fill total cross section into array to avoid repeating calculation ----*/
	const double nuEne_min    = GetEnergyMin();
	const double nuEne_max    = GetEnergyMax();
  const int nuEneNBins      = GetEnergyNBins();
  const double nuEneBinSize = GetEnergyBinWidth();
  const int tNBins          = GetTimeNBins();
  const double tBinSize     = GetTimeBinWidth();
  const double tStart       = GetTimeMin();
  const double tEnd         = GetTimeMax();
  std::vector<double> totcrsIBD(nuEneNBins, 0.); // nu_energy -> total-xsec
  std::vector<double> totcrsNue(nuEneNBins, 0.); // nu_energy -> total-xsec
  std::vector<double> totcrsNueb(nuEneNBins, 0.); // nu_energy -> total-xsec
  std::vector<double> totcrsNux(nuEneNBins, 0.); // nu_energy -> total-xsec
  std::vector<double> totcrsNuxb(nuEneNBins, 0.); // nu_energy -> total-xsec
	std::vector<double> Ocrse0[16][7]; // (rctn==0, ix_state==0) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> Ocrse1[16][7]; // (rctn==0, ix_state==1) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> Ocrse2[16][7]; // (rctn==0, ix_state==2) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> Ocrse3[16][7]; // (rctn==0, ix_state==3) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> Ocrse4[16][7]; // (rctn==0, ix_state==4) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> Ocrsp0[16][7]; // (rctn==1, ix_state==0) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> Ocrsp1[16][7]; // (rctn==1, ix_state==1) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> Ocrsp2[16][7]; // (rctn==1, ix_state==2) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> Ocrsp3[16][7]; // (rctn==1, ix_state==3) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> Ocrsp4[16][7]; // (rctn==1, ix_state==4) [ex_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> OcrseSub[5][32]; // (rctn==0) [ix_sate][channel] -> nu_energy -> total-xsec
	std::vector<double> OcrspSub[5][32]; // (rctn==1) [ix_sate][channel] -> nu_energy -> total-xsec
  std::vector<double> OcrsNC[2][14]; // [rctn][ex_state] -> nu_energy -> total-xsec

	/*-----determine SN direction-----*/
  {
    float sdir[3], ra, dec;
    sn_sundir_( m_sn_date, m_sn_time, sdir, & ra, & dec);
    m_sn_dir[0] = sdir[0];
    m_sn_dir[1] = sdir[1];
    m_sn_dir[2] = sdir[2];
  }

  const int flag_event = GetFlagFillEvent();
  const SKSNSimXSecNuElastic::FLAGETHR flag_elastic_thr = flag_event!=0? SKSNSimXSecNuElastic::ETHROFF : SKSNSimXSecNuElastic::ETHRON;
  /*----------------------------------
   * Build up cross section table for each reaction
   *---------------------------------*/
	std::cout << "calculate cross section and fill to array" << std::endl;
  for(int i_nu_ene =0; i_nu_ene < nuEneNBins; i_nu_ene++) {
    const double nu_energy = nuEne_min + ( double(i_nu_ene) + 0.5 ) * nuEneBinSize;
    double crsOx = 0.;
#ifdef DEBUG
    if(i_nu_ene % 10 == 0) std::cout << "Neutrino Energy  " << nu_energy << std::endl;
#endif
    /*----- inverse beta decay -----*/
    constexpr double eEneThr = 5.0;
    if(nu_energy > eEneThr + DeltaM) totcrsIBD[i_nu_ene] = xsecibd.GetCrosssection(nu_energy);
    //SKSNSimTools::DumpDebugMessage(Form("IBD XSEC %.5g MeV -> %.5g cm2", nu_energy, totcrsIBD[i_nu_ene]));

    /*----- electron elastic -----*/

    totcrsNue[i_nu_ene]  = xsecnuela.GetCrosssection(nu_energy,   PDG_ELECTRON_NEUTRINO, flag_elastic_thr);
    totcrsNueb[i_nu_ene] = xsecnuela.GetCrosssection(nu_energy, - PDG_ELECTRON_NEUTRINO, flag_elastic_thr);
    totcrsNux[i_nu_ene]  = xsecnuela.GetCrosssection(nu_energy,   PDG_MUON_NEUTRINO,     flag_elastic_thr); // Nux (here choosing MuNu but this handles as nu_x
    totcrsNuxb[i_nu_ene] = xsecnuela.GetCrosssection(nu_energy, - PDG_MUON_NEUTRINO,     flag_elastic_thr);

		//std::cout << "start calculate cc reaction crosssection" << std::endl; //nakanisi
    /*----- charged current with oxygen -----*/
    // rcn (reaction): nubar or nu
    for(int rcn=0;rcn<2;rcn++){
      for(int state=0;state<5;state++){
        if(state==0){
          for(int ex_state=0;ex_state<3;ex_state++){
            for(int ch=0;ch<7;ch++){
              if(rcn==0){
                //electron neutrino
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrse0[ex_state][ch].push_back(crsOx);
              }
              else if(rcn==1){
                //anti electron neutrino
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrsp0[ex_state][ch].push_back(crsOx);
              }
            }
          }
        }
        else if(state==1){
          for(int ex_state=0;ex_state<15;ex_state++){
            for(int ch=0;ch<7;ch++){
              if(rcn==0){
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrse1[ex_state][ch].push_back(crsOx);
              }
              else if(rcn==1){
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrsp1[ex_state][ch].push_back(crsOx);
              }
            }
          }
        }
        else if(state==2){
          for(int ex_state=0;ex_state<8;ex_state++){
            for(int ch=0;ch<7;ch++){
              if(rcn==0){
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrse2[ex_state][ch].push_back(crsOx);
              }
              else if(rcn==1){
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrsp2[ex_state][ch].push_back(crsOx);
              }
            }
          }
        }
        else if(state==3){
          for(int ex_state=0;ex_state<1;ex_state++){
            for(int ch=0;ch<7;ch++){
              if(rcn==0){
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrse3[ex_state][ch].push_back(crsOx);
                //if(ch==0)std::cout << nu_energy << " " << state << " " << ex_state << " " << crsOx << " " << std::endl; //nakanisi
              }
              else if(rcn==1){
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrsp3[ex_state][ch].push_back(crsOx);
              }
            }
          }
        }
        else if(state==4){
          for(int ex_state=0;ex_state<16;ex_state++){
            for(int ch=0;ch<7;ch++){
              if(rcn==0){
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrse4[ex_state][ch].push_back(crsOx);
              }
              else if(rcn==1){
                crsOx = xsecnuoxygen.GetCrosssection(nu_energy, {rcn, state, ex_state, ch});
                Ocrsp4[ex_state][ch].push_back(crsOx);
              }
            }
          }
        }
      }
    }
    //calculate cross section of sub channel and excited state
    for(int rcn=0;rcn<2;rcn++){
      for(int state=0;state<5;state++){
        for(int ch=0;ch<32;ch++){
          if(rcn==0){
            crsOx = xsecnuoxygensub.GetCrosssection(nu_energy, {rcn, state, ch});
            OcrseSub[state][ch].push_back(crsOx);
            //std::cout << rcn << " " << state << " " << ch << " " << nu_energy << " " << crsOx << std::endl; //nakanisi
          }
          if(rcn==1){
            crsOx = xsecnuoxygensub.GetCrosssection(nu_energy, {rcn, state, ch});
            OcrspSub[state][ch].push_back(crsOx);
          }
        }
      }
    }
    //calculate cross section of nc reaction
    for(int rcn=0;rcn<2;rcn++){
      switch(rcn){
        case 0: //for p + 15N
          for(int ex_state=0;ex_state<xsecnuoxygennc.GetNumEx(rcn);ex_state++){
            crsOx = xsecnuoxygennc.GetCrosssection(nu_energy, {rcn, ex_state});
            OcrsNC[rcn][ex_state].push_back(crsOx);
          }
          break;

        case 1: // for n + 15O
          for(int ex_state=0;ex_state<xsecnuoxygennc.GetNumEx(rcn);ex_state++){
            crsOx = xsecnuoxygennc.GetCrosssection(nu_energy, {rcn, ex_state});
            //crsOx_nc = ocrs_nc -> CsNuOxyNCNue(rcn, nu_energy);
            OcrsNC[rcn][ex_state].push_back(crsOx);
          }
          break;
      }
    }
  }

  const SKSNSIMENUM::NEUTRINOOSCILLATION nuosctype = GetGeneratorNuOscType();
  const auto &osctuple = SKSNSimPhysConst::NuOscProbCollection.at(nuosctype);
  const double oscnue1 = std::get<0>(osctuple);
  const double oscnue2 = std::get<1>(osctuple);
  const double oscneb1 = std::get<2>(osctuple);
  const double oscneb2 = std::get<3>(osctuple);
  const double oscnux1 = std::get<4>(osctuple);
  const double oscnux2 = std::get<5>(osctuple);
  const double oscnxb1 = std::get<6>(osctuple);
  const double oscnxb2 = std::get<7>(osctuple);
  const double RatioTo10kpc = GetSNDistanceRatioTo10kpc();

  // tempolary buffer 
  double rate = 0.0;

  //expected total number of events
  double totNuebarp = 0.;
  double totNueElastic = 0., totNuebarElastic = 0., totNuxElastic = 0., totNuxbarElastic = 0.;
  double totNueO = 0., totNuebarO = 0.;
  double totNueOsub = 0., totNuebarOsub = 0.;
  double totNcNuep = 0., totNcNuebarp = 0., totNcNuxp = 0., totNcNuxbarp = 0., totNcNuen = 0., totNcNuebarn = 0., totNcNuxn = 0., totNcNuxbarn = 0.;
  std::vector<double> totNcNuepCh(8,0.);
  std::vector<double> totNcNuebarpCh(8,0.);
  std::vector<double> totNcNuxpCh(8,0.);
  std::vector<double> totNcNuxbarpCh(8,0.);
  std::vector<double> totNcNuenCh(4,0.);
  std::vector<double> totNcNuebarnCh(4,0.);
  std::vector<double> totNcNuxnCh(4,0.);
  std::vector<double> totNcNuxbarnCh(4,0.);

	/*---- loop ----*/
  std::cout << "start loop in Process" << std::endl; //nakanisi
  double time;
  double nuEne;
  for(Int_t i_time =0; i_time < tNBins; i_time++) {

    time = tStart + (double(i_time)+0.5)*tBinSize; //center value of each bin[s]
    int itime_sn = int(time);

    if(itime_sn > (int)(tEnd * 1000.)){
      exit(0);
    }

    for(int i_nu_ene =0; i_nu_ene < nuEneNBins; i_nu_ene++) {

      const double nu_energy = nuEne_min + ( double(i_nu_ene) + 0.5 ) * nuEneBinSize;

      const double nspcne  = flux.GetFlux(nu_energy, time, SKSNSimFluxModel::FLUXNUE); //Nue
      const double nspcneb = flux.GetFlux(nu_energy, time, SKSNSimFluxModel::FLUXNUEB); //Nuebar
      const double nspcnx  = flux.GetFlux(nu_energy, time, SKSNSimFluxModel::FLUXNUX); //Nux or Nuexbar

      /*----- inverse beta decay -----*/

      rate = Const_p * (oscneb1*nspcneb + oscneb2*nspcnx) * totcrsIBD[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
			//SKSNSimTools::DumpDebugMessage(Form("IBD rate is: time %.2g Enu %.5g TotCRSIBD %.5g Nspcneb x oscneb1 %.5g Nspcnx x oscneb2 %.5g nuEneBinSize %.2g tBinSize %.2g RatioTo10kpc %.2g -> rate %.5g", time, nu_energy, totcrsIBD[i_nu_ene], oscneb1*nspcneb, oscneb2*nspcnx, nuEneBinSize, tBinSize, RatioTo10kpc, rate));
      totNuebarp += rate;
      if(flag_event == 1) {
        auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, 0 /*nReact*/, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
        evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
      }


      /*----- electron elastic -----*/

      rate = Const_e * (oscnue1*nspcne + oscnue2*nspcnx) * totcrsNue[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
      totNueElastic += rate;
      if(flag_event == 1) {
        auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, 1 /*nReact*/, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
        evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
      }

      rate = Const_e * (oscneb1*nspcneb + oscneb2*nspcnx) * totcrsNueb[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
      totNuebarElastic += rate;
      if(flag_event == 1) {
        auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, 2 /*nReact*/, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
        evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
      }

      rate = Const_e * (oscnux1*nspcnx + oscnux2*nspcne) * totcrsNux[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
      totNuxElastic += rate;
      if(flag_event == 1) {
        auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, 3 /*nReact*/, PDG_MUON_NEUTRINO /*nuType*/, rate);
        evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
      }

      rate = Const_e * (oscnxb1*nspcnx + oscnxb2*nspcneb) * totcrsNuxb[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
      totNuxbarElastic += rate;
      if(flag_event == 1) {
        auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, 4 /*nReact*/, - PDG_MUON_NEUTRINO /*nuType*/, rate);
        evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
      }

      /*----- charged current with oxygen -----*/
      for(int ex_energy=0;ex_energy<5;ex_energy++){
        int rcn = 0;
        if(ex_energy==0){
          for(int ex_state=0;ex_state<3;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrse0[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscnue1*nspcne + oscnue2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              //std::cout << time << " " << nu_energy << " " << ex_state << " " << ch << " " << rate << std::endl; //nakanisi
              totNueO += rate;
              if(flag_event == 1){
                //sReact = to_string(rcn) + to_string(ex_energy) + to_string(ex_state) + to_string(ch);
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                //if(ex_state==2 && ch==6)std::cout << "nReact" << " " << nReact << " " << rcn << " " << ex_energy << " " << ex_state << " " << ch << std::endl; //nakanisi
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /* nuType */, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
        if(ex_energy==1){
          //std::cout << "nue ex_energy=1" << std::endl; //nakanisi
          for(int ex_state=0;ex_state<15;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrse1[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscnue1*nspcne + oscnue2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              totNueO += rate;
              if(flag_event == 1){
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
        if(ex_energy==2){
          //std::cout << "nue ex_energy =2" << std::endl; //nakanisi
          for(int ex_state=0;ex_state<8;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrse2[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscnue1*nspcne + oscnue2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              totNueO += rate;
              if(flag_event == 1){
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
        if(ex_energy==3){
          //std::cout << "nue ex_energy = 3" << std::endl; //nakanisi
          for(int ex_state=0;ex_state<1;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrse3[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscnue1*nspcne + oscnue2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              totNueO += rate;
              if(flag_event == 1){
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
        if(ex_energy==4){
          //std::cout << "nue ex_energy = 4" << std::endl; //nakanisi
          for(int ex_state=0;ex_state<16;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrse4[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscnue1*nspcne + oscnue2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              totNueO += rate;
              if(flag_event == 1){
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
      }
      //std::cout << time << " " << nu_energy << "ex_energy loop end" << std::endl; //nakanisi

      for(int ex_energy=0;ex_energy<5;ex_energy++){
        int rcn = 1;
        if(ex_energy==0){
          for(int ex_state=0;ex_state<3;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrsp0[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscneb1+nspcneb + oscneb2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              //std::cout << time << " " << nu_energy << " " << ex_state << " " << ch << " " << rate << std::endl; //nakanisi
              totNuebarO += rate;
              if(flag_event == 1){
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                //std::cout << "MakeEvent" << std::endl;
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
                //std::cout << "end MakeEvent" << std::endl;
              }
            }
          }
        }
        if(ex_energy==1){
          for(int ex_state=0;ex_state<15;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrsp1[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscneb1*nspcneb + oscneb2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              totNuebarO += rate;
              if(flag_event == 1){
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
        if(ex_energy==2){
          for(int ex_state=0;ex_state<8;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrsp2[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscneb1*nspcneb + oscneb2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              totNuebarO += rate;
              if(flag_event == 1){
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
        if(ex_energy==3){
          for(int ex_state=0;ex_state<1;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrsp3[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscneb1*nspcneb + oscneb2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              totNuebarO += rate;
              if(flag_event == 1){
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
        if(ex_energy==4){
          for(int ex_state=0;ex_state<16;ex_state++){
            for(int ch=0;ch<7;ch++){
              double crsOx = Ocrsp4[ex_state][ch].at(i_nu_ene);
              rate = Const_o * (oscneb1*nspcneb + oscneb2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
              totNuebarO += rate;
              if(flag_event == 1){
                const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + (ex_state+1)*10 + (ch+1);
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
      }

      //std::cout << "start sub reaction culculation" << std::endl; //nakanisi
      //nue + O sub reaction
      for(int ex_energy=0;ex_energy<5;ex_energy++){
        int rcn = 0;
        for(int ch=0;ch<32;ch++){
          //std::cout << "OcrseSub" << " " << ex_energy << " " << ch << std::endl; //nakanisi
          double crsOx = OcrseSub[ex_energy][ch].at(i_nu_ene);
          //std::cout << "end OcrseSub" << std::endl; //nakanisi
          rate = Const_o * (oscnue1*nspcne + oscnue2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
          //std::cout << "sub" << " " << time << " " << nu_energy << " " << crsOx << " " << rate << std::endl; //nakanisi
          totNueOsub += rate;
          //if(crsOx != 0.) SKSNSimTools::DumpDebugMessage(Form("NuOxy rate is: time %.2g Enu %.5g Ocrs %.5g Nspcne x oscne1 %.5g Nspcnx x oscne2 %.5g nuEneBinSize %.2g tBinSize %.2g RatioTo10kpc %.2g -> rate %.5g totNueOsub %.5g", time, nu_energy, crsOx, oscnue1*nspcne, oscnue2*nspcnx, nuEneBinSize, tBinSize, RatioTo10kpc, rate, totNueOsub));
          if(flag_event == 1){
            const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 + 3*100 + 9;
            //std::cout << "MakeEvent" << std::endl; //nakanisi
            if(nu_energy > 15.4){
              auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
              evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
            }
            //std::cout << "end MakeEvent" << std::endl; //nakanisi
          }
        }
      }
      //nue_bar + O sub raction
      for(int ex_energy=0;ex_energy<5;ex_energy++){
        int rcn = 1;
        for(int ch=0;ch<32;ch++){
          double crsOx = OcrspSub[ex_energy][ch].at(i_nu_ene);
          rate = Const_o * (oscneb1*nspcneb + oscneb2*nspcnx) * crsOx * nuEneBinSize * tBinSize * RatioTo10kpc;
          //std::cout << "sub" << " " << time << " " << nu_energy << " " << crsOx << " " << " " << Const_o << " " << oscneb1 << " " << nspcneb << " " << oscneb2 << " " << nspcnx << " " << nuEneBinSize << " " << tBinSize << " " << RatioTo10kpc << " " << rate << std::endl; //nakanisi
          totNuebarOsub += rate;
          if(flag_event == 1){
            const int nReact = (rcn+1)*10e4 + (ex_energy+1)*10e3 +3*100 + 9;
            //if(ch==0 && ex_energy==1)std::cout << "nReact" << " " << nReact << " " << rcn << " " << ex_energy << " " << ch << std::endl; //nakanisi
            if(nu_energy > 11.4) {
              auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
              evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
            }

          }
        }
      }

      // NC reaction
      auto getNeutrinoType = [](int rcn){
        const static int neutrinoType[4] = { PDG_ELECTRON_NEUTRINO, - PDG_ELECTRON_NEUTRINO, PDG_MUON_NEUTRINO, - PDG_MUON_NEUTRINO};
        return neutrinoType[rcn];
      };

      for(int rcn=0;rcn<4;rcn++){
        for(int excit=0;excit<2;excit++){
          if(excit==0){ //p + 15N reaction
            for(int ex_energy=0;ex_energy<8;ex_energy++){
              double crsOx_nc = OcrsNC[0][ex_energy].at(i_nu_ene);
              //double crsOx_nc = OcrsNC[0].at(i_nu_ene);
              if(rcn==0){
                rate = Const_o * (oscnue1*nspcne + oscnue2*nspcnx) * crsOx_nc * nuEneBinSize * tBinSize * RatioTo10kpc;
                totNcNuep += rate;
                totNcNuepCh[ex_energy] += rate;
              }
              else if(rcn==1){
                rate = Const_o * (oscneb1*nspcneb + oscneb2*nspcnx) * crsOx_nc * nuEneBinSize * tBinSize * RatioTo10kpc;
                totNcNuebarp += rate;
                totNcNuebarpCh[ex_energy] += rate;
              }
              else if(rcn==2){
                rate = Const_o * (oscnux1*nspcnx + oscnux2*nspcne) * crsOx_nc * nuEneBinSize * tBinSize * RatioTo10kpc;
                totNcNuxp += rate;
                totNcNuxpCh[ex_energy] += rate;
              }
              else if(rcn==3){
                rate = Const_o * (oscnxb1*nspcnx + oscnxb2*nspcneb) * crsOx_nc * nuEneBinSize * tBinSize * RatioTo10kpc;
                totNcNuxbarp += rate;
                totNcNuxbarpCh[ex_energy] += rate;
              }
              //std::cout << "NC rate: " << time << " " << nu_energy << " " << rcn << " " << excit << " " << crsOx_nc << " " << rate << std::endl; // nakanisi
              //totNcNup += rate;
              if(flag_event == 1){
                const int nReact = 3000 + (rcn+1)*100 + (excit+1)*10 + (ex_energy+1);
                //std::cout << "NC event " << nReact << " " << 3000 << " " << rcn << " " << excit << " " << ex_energy << std::endl;
                const int nuType = getNeutrinoType(rcn);
                auto buf = MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, nuType, rate);

                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
          if(excit==1){ //n + 15O reaction
            for(int ex_energy=0;ex_energy<4;ex_energy++){
              double crsOx_nc = OcrsNC[1][ex_energy].at(i_nu_ene);
              //double crsOx_nc = OcrsNC[1].at(i_nu_ene);
              if(rcn==0){
                rate = Const_o * (oscnue1*nspcne + oscnue2*nspcnx) * crsOx_nc * nuEneBinSize * tBinSize * RatioTo10kpc;
                totNcNuen += rate;
                totNcNuenCh[ex_energy] += rate;
              }
              else if(rcn==1){
                rate = Const_o * (oscneb1*nspcneb + oscneb2*nspcnx) * crsOx_nc * nuEneBinSize * tBinSize * RatioTo10kpc;
                totNcNuebarn += rate;
                totNcNuebarnCh[ex_energy] += rate;
              }
              else if(rcn==2){
                rate = Const_o * (oscnux1*nspcnx + oscnux2*nspcne) * crsOx_nc * nuEneBinSize * tBinSize * RatioTo10kpc;
                totNcNuxn += rate;
                totNcNuxnCh[ex_energy] += rate;
              }
              else if(rcn==3){
                rate = Const_o * (oscnxb1*nspcnx + oscnxb2*nspcneb) * crsOx_nc * nuEneBinSize * tBinSize * RatioTo10kpc;
                totNcNuxbarn += rate;
                totNcNuxbarnCh[ex_energy] += rate;
              }
              //totNcNun += rate;
              if(flag_event == 1){
                const int nReact = 3000 + (rcn+1)*100 + (excit+1)*10 + (ex_energy+1);
                //std::cout << nReact << " rcn " << rcn << " excit " << excit << " ex_energy " << ex_energy << std::endl;
                const int nuType = getNeutrinoType(rcn);
                auto buf =  MakeEvent(nuEneBinSize, tBinSize, time, nu_energy, nReact, nuType, rate);
                evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
              }
            }
          }
        }
      }
    }

    //std::cout << time << " " << totNuebarp << " " << totNueElastic << std::endl;

  }
  std::cout << "end loop process" << std::endl; //nakanisi

  const double totalNumOfEvts = ( totNuebarp + totNueElastic
      + totNuebarElastic + totNuxElastic + totNuxbarElastic
      + totNueO + totNuebarO 
      + totNueOsub + totNuebarOsub
      + totNcNuep + totNcNuen + totNcNuebarp + totNcNuebarn + totNcNuxp + totNcNuxn + totNcNuxbarp + totNcNuxbarn
      );


  fprintf( stdout, "------------------------------------\n" );
  fprintf( stdout, "total expected number of events %e\n", totalNumOfEvts );
  fprintf( stdout, "   nuebar + p = %e\n", totNuebarp );
  fprintf( stdout, "   nue + e = %e\n", totNueElastic );
  fprintf( stdout, "   nuebar + e = %e\n", totNuebarElastic );
  fprintf( stdout, "   nux + e = %e\n", totNuxElastic );
  fprintf( stdout, "   nuxbar + e = %e\n", totNuxbarElastic );
  fprintf( stdout, "   nue + O (CC) = %e\n", totNueO+totNueOsub );
  fprintf( stdout, "   nuebar + O (CC) = %e\n", totNuebarO+totNuebarOsub );
  fprintf( stdout, "   nue + O (NC: p+15N) = %e\n", totNcNuep );
  fprintf( stdout, "   (NC: p+15N) %e, %e, %e, %e, %e, %e, %e, %e\n", totNcNuepCh[0], totNcNuepCh[1], totNcNuepCh[2], totNcNuepCh[3], totNcNuepCh[4], totNcNuepCh[5], totNcNuepCh[6], totNcNuepCh[7] );
  fprintf( stdout, "   nue + O (NC: n+15O) = %e\n", totNcNuen );
  fprintf( stdout, "   (NC: n+15O) %e, %e, %e, %e\n", totNcNuenCh[0], totNcNuenCh[1], totNcNuenCh[2], totNcNuenCh[3] );
  fprintf( stdout, "   nuebar + O (NC: p+15N) = %e\n", totNcNuebarp);
  fprintf( stdout, "   (NC: p+15N) %e, %e, %e, %e, %e, %e, %e, %e\n", totNcNuebarpCh[0], totNcNuebarpCh[1], totNcNuebarpCh[2], totNcNuebarpCh[3], totNcNuebarpCh[4], totNcNuebarpCh[5], totNcNuebarpCh[6], totNcNuebarpCh[7] );
  fprintf( stdout, "   nuebar + O (NC: n+15O) = %e\n", totNcNuebarn);
  fprintf( stdout, "   (NC: n+15O) %e, %e, %e, %e\n", totNcNuebarnCh[0], totNcNuebarnCh[1], totNcNuebarnCh[2], totNcNuebarnCh[3] );
  fprintf( stdout, "   nux + O (NC: p+15N) = %e\n", totNcNuxp);
  fprintf( stdout, "   (NC: p+15N) %e, %e, %e, %e, %e, %e, %e, %e\n", totNcNuxpCh[0], totNcNuxpCh[1], totNcNuxpCh[2], totNcNuxpCh[3], totNcNuxpCh[4], totNcNuxpCh[5], totNcNuxpCh[6], totNcNuxpCh[7] );
  fprintf( stdout, "   nux + O (NC: n+15O) = %e\n", totNcNuxn);
  fprintf( stdout, "   (NC: n+15O) %e, %e, %e, %e\n", totNcNuxnCh[0], totNcNuxnCh[1], totNcNuxnCh[2], totNcNuxnCh[3] );
  fprintf( stdout, "   nuxbar + O (NC: p+15N) = %e\n", totNcNuxbarp);
  fprintf( stdout, "   (NC: p+15N) %e, %e, %e, %e, %e, %e, %e, %e\n", totNcNuxbarpCh[0], totNcNuxbarpCh[1], totNcNuxbarpCh[2], totNcNuxbarpCh[3], totNcNuxbarpCh[4], totNcNuxbarpCh[5], totNcNuxbarpCh[6], totNcNuxbarpCh[7] );
  fprintf( stdout, "   nuxbar + O (NC: n+15O) = %e\n", totNcNuxbarn);
  fprintf( stdout, "   (NC: n+15O) %e, %e, %e, %e\n", totNcNuxbarnCh[0], totNcNuxbarnCh[1], totNcNuxbarnCh[2], totNcNuxbarnCh[3] );
  fprintf( stdout, "------------------------------------\n" );

  std::cout << "end calculation of each expected event number" << std::endl; //nakanisi

  std::cout << "FillEvent start    ( " << evt_buffer.size()  << " evt)" << std::endl;
  if(flag_event == 1) FillEvent(evt_buffer);
  std::cout << "FillEvent finished ( " << evt_buffer.size()  << " evt)" << std::endl;

  return evt_buffer;
}


std::vector<SKSNSimSNEventVector> SKSNSimVectorSNGenerator::MakeEvent(const double nuEneBinSize, const double tBinSize, const double time, const double nu_energy, const int nReact, const int nuType, const double rate){
  // SKSNSimTools::DumpDebugMessage(Form(" MakeEvent time %.2g nuEne %.2g nReact %d nuType %d rate %.2g", time , nu_energy, nReact, nuType, rate));
  std::vector<SKSNSimSNEventVector> buffer;

  //double totcrsIBD[nuEneNBins] = {0.};
  double dRandTotEvts = randomgenerator->Poisson(rate);
  //if(time<0.005)std::cout << time << " " << nu_energy << " " << nReact << " " << nuType << " " << rate << std::endl; //nakanisi

  auto getRandomReal = [](double s, double e , TRandom& rng){ return rng.Uniform(s,e);};

  auto determinePosition = std::bind([](TRandom &rng, SKSNSIMENUM::TANKVOLUME t)
  {
    double rPositionRange = RINTK;
    double hPositionRange = ZPINTK;
    const double FVCUT = 200.;
    switch (t)
    {
      case SKSNSIMENUM::TANKVOLUME::kIDFV: //Fiducial volume
        rPositionRange = RINTK - FVCUT;
        hPositionRange = ZPINTK - FVCUT;
        break;
      case SKSNSIMENUM::TANKVOLUME::kIDFULL: //entire ID volume
        rPositionRange = RINTK;
        hPositionRange = ZPINTK;
        break;
      case SKSNSIMENUM::TANKVOLUME::kTANKFULL: //entire detector volume (including OD)
        rPositionRange = DITKTK;
        hPositionRange = ZPTKTK;
        break;
      default: //entire ID volume
        rPositionRange = RINTK;
        hPositionRange = ZPINTK;
    }
    //random inside the full tank (32.5kton)
    const double r2 = rng.Uniform(1.) * rPositionRange*rPositionRange ;
    const double r = std::sqrt( r2 );
    const double phi = rng.Uniform( 2. * M_PI);
    const double x = r * std::cos( phi );
    const double y = r * std::sin( phi );
    const double z = -hPositionRange + rng.Uniform( 2.*hPositionRange);

    return UtilVector3<double>(x,y,z);
  }, std::placeholders::_1, m_generator_volume);

  if(dRandTotEvts > 0){
    //std::cout << "MakeEvent " << time << " " << nu_energy << " " << nReact << " " << nuType << " " << rate << " " << dRandTotEvts << std::endl; //nakanisi
    for(int i=0; i<dRandTotEvts; i++){
      double ene_s = nu_energy - nuEneBinSize/2., ene_e = nu_energy + nuEneBinSize/2.;
      double nuEne = getRandomReal( ene_s, ene_e, *randomgenerator );
      double time_s = time - tBinSize/2., time_e = time + tBinSize/2.;
      double tReact = randomgenerator->Uniform(time_s, time_e); //getRandomReal( time_s, time_e , randomgenerator );
      //std::cout << tReact << " " << nuEne << " " << nReact << " " << nuType << " " << rate << std::endl; //nakanisi

      double ver_x = 9999., ver_y = 9999., ver_z = 9999.;
      //auto xyz = determinePosition(mInnerID, ver_x, ver_y, ver_z );
      auto xyz = determinePosition(*randomgenerator);

      //SNEvtInfo evtInfo;
      //evtInfo.iEvt = dRandTotEvts;
      SKSNSimSNEventVector evtInfo;

      const double rvtx [3] = {xyz.x, xyz.y, xyz.z};
      evtInfo.SetSNEvtInfo(nReact, tReact, nuType, nuEne, m_sn_dir, rvtx);
      evtInfo.SetRunnum(GetRUNNUM());
      evtInfo.SetSubRunnum(GetSubRUNNUM());

      buffer.push_back( evtInfo );
    }
  }
#ifdef DEBUG
  ///std::cout << "[IZU] MakeEvent returning " << buffer.size() << " events" << std::endl;
#endif
  return buffer;
}

void SKSNSimVectorSNGenerator::FillEvent(std::vector<SKSNSimSNEventVector> &evt_buffer)
{

	/*---- Time sorting ----*/
  std::cout << "start time sorting before loop of FillEvent" << std::endl;
  std::sort( evt_buffer.begin(), evt_buffer.end());

  int totGenNuebarp=0;
  int totGenNueElastic=0, totGenNuebarElastic=0, totGenNuxElastic=0, totGenNuxbarElastic=0;
  int totGenNueO=0, totGenNuebarO=0;
  int totGenNueOsub=0, totGenNuebarOsub=0;
  int totGenNcNuep=0, totGenNcNuebarp=0, totGenNcNuxp=0, totGenNcNuxbarp=0, totGenNcNuen=0, totGenNcNuebarn=0, totGenNcNuxn=0, totGenNcNuxbarn=0;
  std::vector<int> totGenNcNuepCh(8, 0);
  std::vector<int> totGenNcNuebarpCh(8,0);
  std::vector<int> totGenNcNuxpCh(8,0);
  std::vector<int> totGenNcNuxbarpCh(8,0);
  std::vector<int> totGenNcNuenCh(4,0);
  std::vector<int> totGenNcNuebarnCh(4,0);
  std::vector<int> totGenNcNuxnCh(4,0);
  std::vector<int> totGenNcNuxbarnCh(4,0);

  int iSkip = 0;

  std::cout << "start event loop in FillEvent" << std::endl; //nakanisi
  for( uint iEvt = 0; iEvt < evt_buffer.size(); iEvt++ ){

    iSkip = 0;

    SKSNSimSNEventVector & p = evt_buffer[iEvt];

    // fill SNEvtInfo (see $SKOFL_ROOT/include/lowe/snevtinfo.h )

    p.SetSNEvtInfoIEvt(iEvt);

    // MCVERTEX (see $SKOFL_ROOT/inc/vcvrtx.h )

    p.AddVertex(
        p.GetSNEvtInfoRVtx(0),
        p.GetSNEvtInfoRVtx(1),
        p.GetSNEvtInfoRVtx(2),
        1 /* iflvvc */,
        0 /* iparvc */,
        0. /* timvvc */);
    // impossible store here because it is float and no enough precision for SN time, instead of this, fill it into sngen->rTime above


    //std::cout << iEvt << " t=" << p.rTime << " " << p.rType << " " << p.nuType << " E=" << p.nuEne << " x=" << p.rVtx[0] << " y=" << p.rVtx[1] << " z=" << p.rVtx[2] << std::endl;

    // Calculate neutrino interaction vector and save into MCVECT
    determineKinematics( xsecmodels, *randomgenerator, p, m_sn_dir);

    const auto rType = p.GetSNEvtInfoRType();
    if(iSkip == 0) {

      if(rType == 0) totGenNuebarp++;
      else if(rType == 1) totGenNueElastic++;
      else if(rType == 2) totGenNuebarElastic++;
      else if(rType == 3) totGenNuxElastic++;
      else if(rType == 4) totGenNuxbarElastic++;
      else if(rType>1000 && rType<10000){ // nc reaction
        const int Reaction = rType/1000;
        const int Excit_pre = rType/100;
        const int Excit = ((rType - (Reaction)*1000)/100) - 1;
        const int ch_pre = rType/10;
        const int channel = (rType - ch_pre*10) - 1;
        const int particle = ((rType - Excit_pre*100)/10) - 1;
        //std::cout << "NC reaction " << p.rType << " Reaction " << Reaction << " Excit " << Excit << " particle " << particle << " channel " << channel << std::endl; //nakanisi
        if(Excit==0 && particle==0)totGenNcNuep++;
        else if(Excit==1 && particle==0)totGenNcNuebarp++;
        else if(Excit==2 && particle==0)totGenNcNuxp++;
        else if(Excit==3 && particle==0)totGenNcNuxbarp++;
        else if(Excit==0 && particle==1)totGenNcNuen++;
        else if(Excit==1 && particle==1)totGenNcNuebarn++;
        else if(Excit==2 && particle==1)totGenNcNuxn++;
        else if(Excit==3 && particle==1)totGenNcNuxbarn++;
      }
      else if(rType>=10000){ // cc reaction
                               //if(Ex_state==30)std::cout << p.rType << " " << "nReact" << " " << nReact << " " << "Reaction" << " " << Reaction << " " << "State_pre" << " " << State_pre << " " << "State" << " " << State << " " << "Ex_state_pre" << " " << Ex_state_pre << " " << "Ex_state" << " " << Ex_state << " " << "channel" << " " << channel << std::endl; //nakanisi
        const int Reaction = rType/10e4 - 1;
        const int State_pre = rType/10e3;
        const int Ex_state_pre = rType/10;
        const int State = ((rType - (Reaction+1)*10e4)/10e3) - 1;
        const int Ex_state = ((rType - State_pre*10e3)/10) - 1;
        const int channel = (rType - Ex_state_pre*10) - 1;
        if(Reaction==0 && Ex_state!=29) totGenNueO++;
        else if(Reaction==1 && Ex_state!=29) totGenNuebarO++;
        else if(Reaction==0 && Ex_state==29) totGenNueOsub++;
        else if(Reaction==1 && Ex_state==29) totGenNuebarOsub++; 
      }
    }
  }

  const int totalNumOfGenEvts = ( totGenNuebarp + totGenNueElastic
      + totGenNuebarElastic + totGenNuxElastic + totGenNuxbarElastic
      + totGenNueO + totGenNuebarO 
      + totGenNueOsub + totGenNuebarOsub
      + totGenNcNuep + totGenNcNuebarp + totGenNcNuxp + totGenNcNuxbarp + totGenNcNuen + totGenNcNuebarn + totGenNcNuxn + totGenNcNuxbarn
      );

  fprintf( stdout, "------------------------------------\n" );
  fprintf( stdout, "total generated number of events %d\n", totalNumOfGenEvts );
  fprintf( stdout, "   nuebar + p = %d\n", totGenNuebarp );
  fprintf( stdout, "   nue + e = %d\n", totGenNueElastic );
  fprintf( stdout, "   nuebar + e = %d\n", totGenNuebarElastic );
  fprintf( stdout, "   nux + e = %d\n", totGenNuxElastic );
  fprintf( stdout, "   nuxbar + e = %d\n", totGenNuxbarElastic );
  fprintf( stdout, "   nue + o = %d\n", totGenNueO+totGenNueOsub );
  fprintf( stdout, "   nuebar + o = %d\n", totGenNuebarO+totGenNuebarOsub );
  fprintf( stdout, "   nue + o (NC:p+15N) = %d\n", totGenNcNuep );
  fprintf( stdout, "   nuebar + o (NC:p+15N) = %d\n", totGenNcNuebarp );
  fprintf( stdout, "   nux + o (NC:p+15N) = %d\n", totGenNcNuxp );
  fprintf( stdout, "   nuxbar + o (NC:p+15N) = %d\n", totGenNcNuxbarp );
  fprintf( stdout, "   nue + o (NC:n+15O) = %d\n", totGenNcNuen );
  fprintf( stdout, "   nuebar + o (NC:n+15O) = %d\n", totGenNcNuebarn );
  fprintf( stdout, "   nux + o (NC:n+15O) = %d\n", totGenNcNuxn );
  fprintf( stdout, "   nuxbar + o (NC:n+15O) = %d\n", totGenNcNuxbarn );
  fprintf( stdout, "------------------------------------\n" );

}

void SKSNSimVectorSNGenerator::determineKinematics( std::map<XSECTYPE, std::shared_ptr<SKSNSimCrosssectionModel>> xsecmodels, TRandom &rng, SKSNSimSNEventVector &ev, const double snDir[])
{
  auto SQ = [](double x){return x*x;};
  const double nuEne = ev.GetSNEvtInfoNuEne();
  const auto pvect = nuEne * ev.GetSNEvtInfoNuDir();

  auto generateNormVect = [](TRandom &rng)
  {
    //random reaction of neutron
    double phi = rng.Uniform(0., 2.*M_PI); 
    double theta = rng.Uniform(0., M_PI); // TODO this should be dCosTheta? This is based on original code
    return UtilVector3<double>(theta, phi).Unit();
  };

  //number of particle emitted on deexcitation with CC reaction
  constexpr int numNtNueO[7] = {0, 1, 0, 2, 0, 0, 0};
  constexpr int numPtNueO[7] = {1, 1, 2, 1, 0, 0, 1};
  constexpr int numNtNuebarO[7] = {0, 1, 0, 2, 1, 0, 0};
  constexpr int numPtNuebarO[7] = {0, 0, 1, 0, 1, 1, 2};
  constexpr int numGmNuebarO[7] = {1, 0, 0, 0, 0, 0, 0};

  const UtilVector3<double> snDir_vec(snDir);
  int iSkip = 0;

  SKSNSimXSecIBDSV       &xsecibd         = dynamic_cast<SKSNSimXSecIBDSV&>(      *xsecmodels[XSECTYPE::mXSECIBD]);
  SKSNSimXSecNuElastic   &xsecnuela       = dynamic_cast<SKSNSimXSecNuElastic&>(  *xsecmodels[XSECTYPE::mXSECELASTIC]);
  SKSNSimXSecNuOxygen    &xsecnuoxygen    = dynamic_cast<SKSNSimXSecNuOxygen&>(   *xsecmodels[XSECTYPE::mXSECOXYGEN]);
  SKSNSimXSecNuOxygenSub &xsecnuoxygensub = dynamic_cast<SKSNSimXSecNuOxygenSub&>(*xsecmodels[XSECTYPE::mXSECOXYGENSUB]);

	double sn_theta = acos( snDir[2] );
	double sn_phi = atan2( snDir[1],  snDir[0] );
  const UtilMatrix3<double> Rmat( cos(sn_theta)*cos(sn_phi), -sin(sn_phi), sin(sn_theta)*cos(sn_phi),
                                  cos(sn_theta)*sin(sn_phi),  cos(sn_phi), sin(sn_theta)*sin(sn_phi),
                                             -sin(sn_theta),           0.,             cos(sn_theta));

  auto nReact = ev.GetSNEvtInfoRType();
  if( nReact == 0 ){ // nuebar + p -> e+ + n
    const auto nuMomentum = pvect;
    determineKinematicsIBD( xsecibd, rng, ev, nuMomentum);

  } else if( nReact == 1 || nReact == 2 || nReact == 3 || nReact == 4 ){ //nu + e Elastic
                                                                         //mc->mcinfo[0] = 85007;
                                                                         // Original neutrino
    int ipvc_tmp = 0;
    if( nReact == 1) ipvc_tmp =  PDG_ELECTRON_NEUTRINO;
    if( nReact == 2) ipvc_tmp = -PDG_ELECTRON_NEUTRINO;
    if( nReact == 3) ipvc_tmp =  PDG_MUON_NEUTRINO;
    if( nReact == 4) ipvc_tmp = -PDG_MUON_NEUTRINO;
    auto mom = nuEne * snDir_vec;
    ev.AddTrack(
        ipvc_tmp, nuEne,
        mom.x, mom.y, mom.z,
        0 /*iorgvc*/,
        1 /*ivtivc*/,
        1 /*ivtfvc*/,
        -1 /*iflgvc*/,
        0 /*icrnvc*/
        );

    // Recoil electron
    double eEne, eTheta, ePhi;
    determineAngleElastic( rng, xsecnuela, nReact, nuEne, eEne, eTheta, ePhi, iSkip);
    double amom = sqrt(SQ( eEne ) - SQ( Me ));

    const UtilVector3<double> eDir = Rmat * UtilVector3<double>(eTheta,ePhi); 
    auto eleMom = amom * eDir;
    ev.AddTrack(
        PDG_ELECTRON, eEne,
        eleMom.x, eleMom.y, eleMom.z,
        1 /*iorgvc*/,
        1 /*ivtivc*/,
        1 /*ivtfvc*/,
        0 /*iflgvc*/,
        1 /*icrnvc*/
        );

    double costh = snDir_vec * eDir; //[0] * eDir[0] + snDir[1] * eDir[1] + snDir[2] * eDir[2];
    //std::cout << "Elastic " << costh << std::endl;

  } else if(nReact>1000 && nReact < 10000){
    // Oxygen NC

    //energy of gamma on deexcitation with NC reaction
    const static double eneGamN[8] = {5.27, 6.33, 7.16, 7.56, 8.32, 8.57, 9.05, 9.76};
    const static double eneGamO[4] = {5.18, 6.18, 6.69, 7.28};

    const int Reaction = nReact/1000;
    const int Excit_pre = nReact/100;
    const int Excit = ((nReact - (Reaction)*1000)/100) - 1;
    const int ch_pre = nReact/10;
    const int channel = (nReact - ch_pre*10) - 1;
    const int particle = ((nReact - Excit_pre*100)/10) - 1;

    const auto nuMomentum = pvect;

    int tmp_ipvc = 0;
    if(Excit == 0)      tmp_ipvc =   PDG_ELECTRON_NEUTRINO /* 12 */;
    else if(Excit == 1) tmp_ipvc = - PDG_ELECTRON_NEUTRINO /* -12 */;
    else if(Excit == 2) tmp_ipvc =   PDG_MUON_NEUTRINO /* 14 */;
    else if(Excit == 3) tmp_ipvc = - PDG_MUON_NEUTRINO /* -14 */;

    ev.AddTrack(
        tmp_ipvc, nuEne,
        nuMomentum.x, nuMomentum.y, nuMomentum.z,
        0 /* iorgvc */,
        1 /* ivtivc */,
        1 /* ivtfvc */,
        -1 /* iflgvc */,
        0 /* icrnvc */
        );

    double pTheta, pPhi, pDir[3];
    if(particle == 0){
      double amom = sqrt(SQ(Mp+0.5) - SQ(Mp));
      const auto protonMomentum = amom * generateNormVect(rng);
      ev.AddTrack(
          PDG_PROTON, Mp + 0.5,
          protonMomentum.x, protonMomentum.y, protonMomentum.z,
          1 /* iorgvc */,
          1 /* ivtivc */,
          1 /* ivtfvc */,
          0 /* iflgvc */,
          1 /* icrnvc */
          );

      // Gamma ray
      const auto gammaMomentum = eneGamN[channel] * generateNormVect(rng);
      ev.AddTrack(
          PDG_GAMMA, eneGamN[channel],
          gammaMomentum.x, gammaMomentum.y, gammaMomentum.z,
          1 /* iorgvc */,
          1 /* ivtivc */,
          1 /* ivtfvc */,
          0 /* iflgvc */,
          1 /* icrnvc */
          );
      std::cout << "NC gamma emission " << PDG_GAMMA << " " << eneGamN[channel] << " " << gammaMomentum.x << " " << gammaMomentum.y << " " << gammaMomentum.z << std::endl; // nakanisi
    }
    else if(particle == 1){
      const double amom = sqrt(SQ(0.5+Mn) - SQ(Mn));
      const auto neutronMomentum = amom * generateNormVect(rng);
      ev.AddTrack(
          PDG_NEUTRON, 0.5 + Mn,
          neutronMomentum.x, neutronMomentum.y, neutronMomentum.z,
          1 /* iorgvc */,
          1 /* ivtivc */,
          1 /* ivtfvc */,
          0 /* iflgvc */,
          1 /* icrnvc */
          );

      // Gamma ray
      const auto gammaMomentum = eneGamO[channel] * generateNormVect(rng);
      ev.AddTrack(
          PDG_GAMMA, eneGamO[channel],
          gammaMomentum.x, gammaMomentum.y, gammaMomentum.z,
          1 /* iorgvc */,
          1 /* ivtivc */,
          1 /* ivtfvc */,
          0 /* iflgvc */,
          1 /* icrnvc */
          );
                         //std::cout << "gamma emission " << i_nucre << " " << mc->ipvc[mc->nvc] << " " << mc->energy[mc->nvc] << " " << x << " " << y << " " << z << std::endl; // nakanisi
      std::cout << "NC gamma emission " << PDG_GAMMA << " " << eneGamO[channel] << " " << gammaMomentum.x << " " << gammaMomentum.y << " " << gammaMomentum.z << std::endl; // nakanisi
    }
  } else if(nReact >= 10000){
    // Oxygen CC
    //mc->mcinfo[0] = 85005;
    const int Reaction = nReact/10e4 - 1;
    const int State_pre = nReact/10e3;
    const int Ex_state_pre = nReact/10;
    const int State = ((nReact - (Reaction+1)*10e4)/10e3) - 1;
    const int Ex_state = ((nReact - State_pre*10e3)/10) - 1;
    const int channel = (nReact - Ex_state_pre*10) - 1;
    //if(Reaction==0 && Ex_state!=29)mc->nvc = 2 + numNtNueO[channel];
    //if(Reaction==1 && Ex_state!=29 && channel!=0)mc->nvc = 2 + numNtNuebarO[channel];
    //if(Reaction==1 && Ex_state!=29 && channel==0)mc->nvc = 2 + numGmNuebarO[channel];
    //
    if(Ex_state==29){
      //std::cout << "nReact" << " " << nReact << " " << "Reaction" << " " << Reaction << " " << "State_pre" << " " << State_pre << " " << "State" << " " << State << " " << "Ex_state_pre" << " " << Ex_state_pre << " " << "Ex_state" << " " << Ex_state << " " << "channel" << " " << channel << " " << mc->nvc << std::endl; //nakanisi
    }
    // Original neutrino
    //if(Ex_state==29)mc->nvc = 2;
    double eEne, eTheta, ePhi;
    determineAngleNueO( rng, xsecnuoxygen, Reaction, State, Ex_state, channel, nuEne, eEne, eTheta, ePhi); // channel = 8 is sub reaction of NueO

    int ipvc_tmp = 0;
    auto nuMom = nuEne * snDir_vec;
    if(Reaction == 0) ipvc_tmp =  PDG_ELECTRON_NEUTRINO;
    if(Reaction == 1) ipvc_tmp = -PDG_ELECTRON_NEUTRINO;
    ev.AddTrack(
        ipvc_tmp, nuEne,
        nuMom.x, nuMom.y, nuMom.z,
        0 /*iorgvc*/,
        1 /*ivtivc*/,
        1 /*ivtfvc*/,
        -1 /*iflgvc*/,
        0 /*icrnvc*/
        );
    ev.AddTrack( // Oxygen
        1000080160, 0.,
        0., 0., 0.,
        0 /*iorgvc*/,
        1 /*ivtivc*/,
        1 /*ivtfvc*/,
        -1 /*iflgvc*/,
        0 /*icrnvc*/
        );

    // electron or positron
    //double eEne, eTheta, ePhi, eDir[3];
    //determineAngleNueO( Reaction, State, Ex_state, channel, nuEne, eEne, eTheta, ePhi); // channel = 8 is sub reaction of NueO
    //if(Ex_state==29)std::cout << "determine eEnergy " << Reaction << " " << State << " " << Ex_state << " " << channel << " " << nuEne << " " << eEne << " " << eTheta << " " << ePhi << std::endl; // nakanisi
    double amom = sqrt(SQ( eEne ) - SQ( Me ));
    const UtilVector3<double> eDir = Rmat * UtilVector3<double>(eTheta, ePhi);

    if(Reaction==0) ipvc_tmp =  PDG_ELECTRON; // electron
    if(Reaction==1) ipvc_tmp = -PDG_ELECTRON; // positron
    auto eleMom = amom * UtilVector3<double>(eDir);
    ev.AddTrack(
        ipvc_tmp, eEne,
        eleMom.x, eleMom.y, eleMom.z,
        1 /*iorgvc*/,
        1 /*ivtivc*/,
        1 /*ivtfvc*/,
        0 /*iflgvc*/,
        1 /*icrnvc*/
        );

    //if(eEne<0.)std::cout << "e-/e+ momentum " << mc->pvc[1][0] << " " << mc->pvc[1][1] << " " << mc->pvc[1][2] << " " << eEne << " " << Me << " " << amom << " " << eDir[0] << " " << eDir[1] << " " << eDir[2] << std::endl; // nakanisi

    double costh = snDir[0] * eDir[0] + snDir[1] * eDir[1] + snDir[2] * eDir[2];

    if(numNtNueO[channel]!=0 || numNtNuebarO[channel]!=0 || numGmNuebarO[channel]!=0){
      if(Reaction==0){
        int i_nucre = 0;
        if(Ex_state!=29){
          for(int i=0;i<numNtNueO[channel];i++){
            // Neutron
            i_nucre++;
            const auto neutronMom = sqrt(SQ(0.5+Mn) - SQ(Mn)) * generateNormVect(rng);
            ev.AddTrack(
                PDG_NEUTRON, 0.5 + Mn,
                neutronMom.x, neutronMom.y, neutronMom.z,
                1 /*iorgvc*/,
                1 /*ivtivc*/,
                1 /*ivtfvc*/,
                0 /*iflgvc*/,
                1 /*icrnvc*/
                );
            //std::cout << "neutron information " << nReact << " " << i_nucre << " " << mc->ipvc[2+i_nucre] << " " << x << " " << y << " " << z << std::endl;
          }
        }
      }
      if(Reaction==1){
        int i_nucre = 0;
        if(Ex_state!=29){
          if(channel!=0){
            for(int i=0;i<numNtNuebarO[channel];i++){
              // Neutron
              i_nucre++;
              auto neutronMom = sqrt(SQ(0.5+Mn)-SQ(Mn)) * generateNormVect(rng);
              ev.AddTrack(
                  PDG_NEUTRON, 0.5 + Mn,
                  neutronMom.x, neutronMom.y, neutronMom.z,
                  1 /*iorgvc*/,
                  1 /*ivtivc*/,
                  1 /*ivtfvc*/,
                  0 /*iflgvc*/,
                  1 /*icrnvc*/
                  );
            }
          }
          if(channel==0){
            // Gamma ray
            i_nucre++;
            constexpr double gammaEne = 12.674;
            const auto gammaMom = gammaEne * generateNormVect(rng);
            ev.AddTrack(
                PDG_GAMMA, gammaEne,
                gammaMom.x, gammaMom.y, gammaMom.z,
                1 /*iorgvc*/,
                1 /*ivtivc*/,
                1 /*ivtfvc*/,
                0 /*iflgvc*/,
                1 /*icrnvc*/
                );
            //std::cout << "gamma emission " << i_nucre << " " << mc->ipvc[mc->nvc] << " " << mc->energy[mc->nvc] << " " << x << " " << y << " " << z << std::endl; // nakanisi
          }
        }
      }
    }
  }
  return;
}

void SKSNSimVectorSNGenerator::determineKinematicsIBD( const SKSNSimXSecIBDSV & xsec, TRandom &rng, SKSNSimSNEventVector &ev, const UtilVector3<double> nuMomentum){
  const UtilVector3<double> nuDir = nuMomentum.Unit();
  const double theta = acos( nuDir[2] );
  const double phi = atan2( nuDir[1],  nuDir[0] );
  const UtilMatrix3<double> Rmat( cos(theta)*cos(phi), -sin(phi), sin(theta)*cos(phi),
      cos(theta)*sin(phi),  cos(phi), sin(theta)*sin(phi),
      -sin(theta),        0.,          cos(theta));

  const double nuEne = nuMomentum.Mag();
  // Original neutrino
  ev.AddTrack(
      - PDG_ELECTRON_NEUTRINO, nuEne,
      nuMomentum.x, nuMomentum.y, nuMomentum.z,
      0 /*iorgvc*/,
      1 /*ivtivc*/,
      1 /*ivtfvc*/,
      -1 /*iflgvc*/,
      0 /*icrnvc*/
      );

  // Original proton
  ev.AddTrack(
      PDG_PROTON, Mp,
      0., 0., 0.,
      0 /*iorgvc*/,
      1 /*ivtivc*/,
      1 /*ivtfvc*/,
      -1 /*iflgvc*/,
      0 /*icrnvc*/
      );

  // Positron
  double eEne, eTheta, ePhi;
  determineAngleNuebarP( rng, xsec, nuEne, eEne, eTheta, ePhi );
  auto SQ = [](double x){return x*x;};
  const double amom = sqrt(SQ( eEne ) - SQ( Me ));
  const UtilVector3<double> eDir = Rmat * UtilVector3<double>( eTheta, ePhi);
  //std::cout << nuEne << " " << eEne << " " << eDir[0] << " " << eDir[1] << " " << eDir[2] << std::endl;

  const auto positronMomentum = amom * eDir.Unit();
  ev.AddTrack(
      - PDG_ELECTRON, eEne,
      positronMomentum.x, positronMomentum.y, positronMomentum.z,
      1 /*iorgvc*/,
      1 /*ivtivc*/,
      1 /*ivtfvc*/,
      0 /*iflgvc*/,
      1 /*icrnvc*/
      );

  // Neutron
  const UtilVector3<double> neutronMomentum = nuMomentum - positronMomentum;
  ev.AddTrack(
      PDG_NEUTRON,
      std::sqrt(neutronMomentum.Mag2() + Mn*Mn),
      neutronMomentum.x, neutronMomentum.y, neutronMomentum.z,
      1 /*iorgvc*/,
      1 /*ivtivc*/,
      1 /*ivtfvc*/,
      0 /*iflgvc*/,
      1 /*icrnvc*/
      );
}

double SKSNSimVectorSNGenerator::FindMaxProb ( const double ene, const SKSNSimCrosssectionModel &xsec){
	//find maximum values, which depends on nuEne
	double maxP = 0.;
  double cost;
  constexpr int costNBins = 1000; // TODO modify configureable bin size
  constexpr double costMin = -1.;
  constexpr double costMax =  1.;
  constexpr double costBinSize = (costMax - costMin) / (double) costNBins;
	for( int iCost = 0; iCost < costNBins; iCost++ ){
		cost = costMin + costBinSize * ( double(iCost) + 0.5 );
		double p = xsec.GetDiffCrosssection(ene, cost).first;

		if( p > maxP ){ maxP = p; }
	}
  return maxP;
}

void SKSNSimVectorSNGenerator::determineAngleNuebarP( TRandom &rng, const SKSNSimXSecIBDSV &xsec, const double nuEne, double & eEne, double & eTheta, double & ePhi )
{

	double nuEnergy = nuEne;
	double cost, eMom, p;
  constexpr double costMin = -1.;
  constexpr double costMax = 1.;

  const double maxP = FindMaxProb(nuEnergy, xsec);

  auto getRandomReal = [](double s, double e , TRandom& rng){ return rng.Uniform(s,e);};

  while( 1 ){
    cost = getRandomReal( costMin, costMax, rng );
    auto ss = xsec.GetDiffCrosssection(nuEnergy, cost);
    eEne = ss.second;
    p = ss.first;

    double x = getRandomReal( 0., maxP, rng );
    if( x < p ){
      eTheta = acos( cost );
      ePhi = getRandomReal( -M_PI,  M_PI, rng );
      break;
    }
  }
	return;
}

void SKSNSimVectorSNGenerator::determineAngleElastic( TRandom &rng, const SKSNSimXSecNuElastic & xsec, const int nReact, const double nuEne, double & eEne, double & eTheta, double & ePhi, int &iSkip )
{

	double nuEnergy = nuEne;
	double cost, eEnergy;
	double p=0., x;

  auto getRandomReal = [](double s, double e , TRandom& rng){ return rng.Uniform(s,e);};

	//we know the maximum prob. happens at cost=1
	cost = 1. - ZERO_PRECISION;
	eEnergy = SKSNSimXSecNuElastic::CalcElectronTotEnergy( nuEnergy, cost );

  //SKSNSimTools::DumpDebugMessage(Form("nReact %d nuEnergy %.5g eEnergy %.5g cost %.5g", nReact, nuEnergy, eEnergy, cost));
	double maxP = 0.;
	if( nReact == 1 ) maxP = sl_nue_dif_rad_( & nuEnergy, & eEnergy);
	if( nReact == 2 ) maxP = sl_neb_dif_rad_( & nuEnergy, & eEnergy);
	if( nReact == 3 ) maxP = sl_num_dif_rad_( & nuEnergy, & eEnergy);
	if( nReact == 4 ) maxP = sl_nmb_dif_rad_( & nuEnergy, & eEnergy);
	maxP *= SKSNSimXSecNuElastic::CalcDeEneDCost( nuEnergy, cost );

	/*
	//take into account the electron energy threshold
	//double costTh = nucrs->calcCosTthElastic( nuEnergy, eEneThr );
	*/

	// avoid too low energy event
	constexpr double eEneThrElastic = 0.51099906 /*Me*/ + ZERO_PRECISION; // TODO Me is different between sl_nue_dif_rad.F and SKSNSimConstant. In order to avoid error in sl_nue_dif_rad.F, I use sl_....F version
	const double costTh = SKSNSimXSecNuElastic::CalcCosThr( nuEnergy, eEneThrElastic );

	if(fabs(costTh) > 1.){ // neutrino energy is too low
	  std::cerr << "Neutrino energy is too low : neutrino energy " << nuEnergy << std::endl;
	  iSkip = 1;
	  return;
	}
  //SKSNSimTools::DumpDebugMessage(Form(" Found maxP determineAngleElastic %.5g", maxP) );

	while( 1 ){
		cost = getRandomReal( costTh, 1., rng);
    eEnergy = SKSNSimXSecNuElastic::CalcElectronTotEnergy( nuEnergy, cost );

		if( nReact == 1 ) p = sl_nue_dif_rad_( & nuEnergy, & eEnergy);
		if( nReact == 2 ) p = sl_neb_dif_rad_( & nuEnergy, & eEnergy);
		if( nReact == 3 ) p = sl_num_dif_rad_( & nuEnergy, & eEnergy);
		if( nReact == 4 ) p = sl_nmb_dif_rad_( & nuEnergy, & eEnergy);
    p *= SKSNSimXSecNuElastic::CalcDeEneDCost( nuEnergy, cost );

		x = getRandomReal( 0., maxP, rng );
		if( x < p ){
			eEne = eEnergy;
			eTheta = acos( cost );
			ePhi = getRandomReal( -M_PI, M_PI, rng );
			break;
		}
	}
  //SKSNSimTools::DumpDebugMessage(" Exit determineAngleElastic " );
  return;
}

void SKSNSimVectorSNGenerator::determineAngleNueO(TRandom &rng, SKSNSimXSecNuOxygen &xsec, const int Reaction, const int State, const int Ex_state, const int channel, const double nuEne, double & eEne, double & eTheta, double & ePhi )
{
	const double nuEnergy = nuEne;
	double cost, p, x, eEnergy;

	// find maximum values, which depends on nuEne
	double maxP = 0.;
	int rcn = 0;
  constexpr int costNBins = 1000;
  constexpr double costMin = -1.;
  constexpr double costMax =  1.;
  constexpr double costBinSize = (costMax - costMin) / (double) costNBins;

  for(int iCost=0;iCost<costNBins;iCost++){
    cost = costMin + costBinSize * (iCost + 0.5);
    p = xsec.OxigFuncAngleRecCC(Reaction, State, Ex_state, channel, nuEne, cost);
    //if(Reaction==0 && State==3 && Ex_state==29 && channel==8) std::cout << "AngleRecCC sub" << " " << p << std::endl;
    eEnergy = xsec.OxigFuncRecEneCC(Reaction, State, Ex_state, channel, nuEne);
    //if(Reaction==0 && State==3 && Ex_state==29 && channel==8) std::cout << "RecEneCC" << " " << eEnergy << std::endl;
    //std::cout << "determineAngleNueO" << " " << Reaction << " " << State << " " << Ex_state << " " << channel << " " << nuEne << " " << cost << " " << p << " " << eEnergy << std::endl; //nakanisi
    if(p>maxP){ maxP = p; } // TODO in this case, maxP is determined by nuEnergy and cos theta both. Is it fine?
    while(1){
      cost = rng.Uniform(costMin, costMax);
      //dir_oxigfunc_( & nuEnergy, & cost, & p, & eEnergy );
      x = rng.Uniform( 0., maxP);
      //std::cout << maxP << " " << p << " " << x << std::endl; //nakanisi
      if(x<p){
        eTheta = acos( cost );
        ePhi = rng.Uniform(-M_PI, M_PI);
        eEne = eEnergy;
        //std::cout << "break" << " " << Reaction << " " << State << " " << Ex_state << " " << channel << " " << eTheta << " " << ePhi << " " << eEnergy << std::endl; //nakanisi
        break;
      }
    }
  }
	return;
}

