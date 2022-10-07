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
  auto determinePosition = std::bind([](TRandom &rng)
  {
    const double rPositionRange = RINTK;
    const double hPositionRange = ZPINTK;
    //random inside the full tank (32.5kton)
    const double r2 = rng.Uniform(1.) * rPositionRange*rPositionRange ;
    const double r = std::sqrt( r2 );
    const double phi = rng.Uniform( 2. * M_PI);
    const double x = r * std::cos( phi );
    const double y = r * std::sin( phi );
    const double z = -hPositionRange + rng.Uniform( 2.*hPositionRange);

    return UtilVector3<double>(x,y,z);
  }, *randomgenerator);
  const auto xyz = determinePosition();

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

SKSNSimVectorSNGenerator::SKSNSimVectorSNGenerator(): m_generator_energy_min(0.0), m_generator_energy_max(300.0){
  xsecmodels[mXSECIBD]       = std::make_unique<SKSNSimXSecIBDSV>();
  xsecmodels[mXSECELASTIC]   = std::make_unique<SKSNSimXSecNuElastic>();
  xsecmodels[mXSECOXYGEN]    = std::make_unique<SKSNSimXSecNuOxygen>();
  xsecmodels[mXSECOXYGENSUB] = std::make_unique<SKSNSimXSecNuOxygenSub>();
}

std::vector<SKSNSimSNEventVector> SKSNSimVectorSNGenerator::GenerateEvents(){
  std::vector<SKSNSimSNEventVector> evt_buffer;
  SKSNSimBinnedFluxModel &flux = dynamic_cast<SKSNSimBinnedFluxModel&>(*fluxmodels[0]); // TODO selectable flux
  if(&flux == NULL) {
    std::cerr << "In GenerateEvents() no appropriate flux model (binned flux)" << std::endl;
    evt_buffer.clear();
    return evt_buffer;
  }

  SKSNSimXSecIBDSV       &xsecibd         = dynamic_cast<SKSNSimXSecIBDSV&>(      *xsecmodels[mXSECIBD]);
  SKSNSimXSecNuElastic   &xsecnuela       = dynamic_cast<SKSNSimXSecNuElastic&>(  *xsecmodels[mXSECELASTIC]);
  SKSNSimXSecNuOxygen    &xsecnuoxygen    = dynamic_cast<SKSNSimXSecNuOxygen&>(   *xsecmodels[mXSECOXYGEN]);
  SKSNSimXSecNuOxygenSub &xsecnuoxygensub = dynamic_cast<SKSNSimXSecNuOxygenSub&>(*xsecmodels[mXSECOXYGENSUB]);

	std::cout << "Prcess of sn_burst side" << std::endl;//nakanisi
	/*---- Fill total cross section into array to avoid repeating calculation ----*/
	const double nuEne_min = GetEnergyMin();
	const double nuEne_max = GetEnergyMax();
  const int nuEneNBins = 3000; //flux.GetNBinsEne();
  const int tNBins = (20.0 - 0.0) / 1.0e-3; //flux.GetNBinsTime();
  const double nuEneBinSize = 300. / 3000.; //flux.GetBinWidthEne(0);
  const double tBinSize = 1.0e-3;//flux.GetBinWidthTime(0);
  constexpr double tStart = 0.;
  constexpr double tEnd = 100000000000.;
  std::vector<double> totcrsIBD(nuEneNBins, 0.);
  std::vector<double> totcrsNue(nuEneNBins, 0.);
  std::vector<double> totcrsNueb(nuEneNBins, 0.);
  std::vector<double> totcrsNux(nuEneNBins, 0.);
  std::vector<double> totcrsNuxb(nuEneNBins, 0.);
	std::vector<double> Ocrse0[16][7];
	std::vector<double> Ocrse1[16][7];
	std::vector<double> Ocrse2[16][7];
	std::vector<double> Ocrse3[16][7];
	std::vector<double> Ocrse4[16][7];
	std::vector<double> Ocrsp0[16][7];
	std::vector<double> Ocrsp1[16][7];
	std::vector<double> Ocrsp2[16][7];
	std::vector<double> Ocrsp3[16][7];
	std::vector<double> Ocrsp4[16][7];
	std::vector<double> OcrseSub[5][32];
	std::vector<double> OcrspSub[5][32];

	/*-----determine SN direction-----*/
	int sn_date[3] = {2011, 3, 23};
  int sn_time[3] = {0, 0, 0,};
	float sdir[3], ra, dec;
	sn_sundir_( sn_date, sn_time, sdir, & ra, & dec);
  sn_dir[0] = sdir[0];
  sn_dir[1] = sdir[1];
  sn_dir[2] = sdir[2];



  constexpr int flag_event = 1;
  constexpr SKSNSimXSecNuElastic::FLAGETHR flag_elastic_thr = flag_event!=0? SKSNSimXSecNuElastic::ETHROFF : SKSNSimXSecNuElastic::ETHRON;
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
  }

  /* tentative constant */
  constexpr double oscnue1 = 1.0;
  constexpr double oscnue2 = 0.0;
  constexpr double oscneb1 = 1.0;
  constexpr double oscneb2 = 0.0;
  constexpr double oscnux1 = 2.0;
  constexpr double oscnux2 = 0.0;
  constexpr double oscnxb1 = 2.0;
  constexpr double oscnxb2 = 0.0;
  constexpr double RatioTo10kpc = 1.0;
  //=== tentative constant 
  double rate = 0.0;
  double totNuebarp = 0.0;
  double totNueElastic = 0.0;
  double totNuebarElastic = 0.0;
  double totNuxElastic = 0.0;
  double totNuxbarElastic = 0.0;
  double totNueO = 0.0;
  double totNuebarO = 0.0;
  double totNueOsub = 0.0;
  double totNuebarOsub = 0.0;
  double totNcNup = 0.0;
  double totNcNun = 0.0;
  double totNcNubarp = 0.0;
  double totNcNubarn = 0.0;

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
        auto buf = MakeEvent(time, nu_energy, 0 /*nReact*/, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
        evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
      }


      /*----- electron elastic -----*/

      rate = Const_e * (oscnue1*nspcne + oscnue2*nspcnx) * totcrsNue[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
      totNueElastic += rate;
      if(flag_event == 1) {
        auto buf = MakeEvent(time, nu_energy, 1 /*nReact*/, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
        evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
      }

      rate = Const_e * (oscneb1*nspcneb + oscneb2*nspcnx) * totcrsNueb[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
      totNuebarElastic += rate;
      if(flag_event == 1) {
        auto buf = MakeEvent(time, nu_energy, 2 /*nReact*/, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
        evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
      }

      rate = Const_e * (oscnux1*nspcnx + oscnux2*nspcne) * totcrsNux[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
      totNuxElastic += rate;
      if(flag_event == 1) {
        auto buf = MakeEvent(time, nu_energy, 3 /*nReact*/, PDG_MUON_NEUTRINO /*nuType*/, rate);
        evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
      }

      rate = Const_e * (oscnxb1*nspcnx + oscnxb2*nspcneb) * totcrsNuxb[i_nu_ene] * nuEneBinSize * tBinSize * RatioTo10kpc;
      totNuxbarElastic += rate;
      if(flag_event == 1) {
        auto buf = MakeEvent(time, nu_energy, 4 /*nReact*/, - PDG_MUON_NEUTRINO /*nuType*/, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /* nuType */, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
                auto buf = MakeEvent(time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
              auto buf = MakeEvent(time, nu_energy, nReact, PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
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
              auto buf = MakeEvent(time, nu_energy, nReact, - PDG_ELECTRON_NEUTRINO /*nuType*/, rate);
              evt_buffer.insert(evt_buffer.end(), buf.begin(), buf.end());
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
      + totNcNup + totNcNun + totNcNubarp + totNcNubarn
      );


  fprintf( stdout, "------------------------------------\n" );
  fprintf( stdout, "total expected number of events %e\n", totalNumOfEvts );
  fprintf( stdout, "   nuebar + p = %e\n", totNuebarp );
  fprintf( stdout, "   nue + e = %e\n", totNueElastic );
  fprintf( stdout, "   nuebar + e = %e\n", totNuebarElastic );
  fprintf( stdout, "   nux + e = %e\n", totNuxElastic );
  fprintf( stdout, "   nuxbar + e = %e\n", totNuxbarElastic );
  fprintf( stdout, "   nue + O = %e\n", totNueO+totNueOsub );
  fprintf( stdout, "   nuebar + O = %e\n", totNuebarO+totNuebarOsub );
  //fprintf( stdout, "   nue + O sub = %e\n", totNueOsub );
  //fprintf( stdout, "   nuebar + O sub = %e\n", totNuebarOsub );
  fprintf( stdout, "------------------------------------\n" );

  std::cout << "end calculation of each expected event number" << std::endl; //nakanisi

#ifdef DEBUG
  fprintf( stdout, "------------------------------------\n" );
  fprintf( stdout, "totreaction:\n" );
  fprintf( stdout, "   nuebar + p = %e\n", totNuebarp );
  fprintf( stdout, "   nue + e = %e\n", totNueElastic );
  fprintf( stdout, "   nuebar + e = %e\n", totNuebarElastic );
  fprintf( stdout, "   nux + e = %e\n", totNuxElastic );
  fprintf( stdout, "   nuxbar + e = %e\n", totNuxbarElastic );
  fprintf( stdout, "   nue + O = %e\n", totNueO );
  fprintf( stdout, "   nuebar + O = %e\n", totNuebarO );
  fprintf( stdout, "   nu + O = nu + p + N%e\n", totNcNup );
  fprintf( stdout, "   nu + O = nu + n + O%e\n", totNcNun );
  fprintf( stdout, "   nubar + O = nubar + p + N%e\n", totNcNubarp );
  fprintf( stdout, "   nubar + O = nubar + n + O%e\n", totNcNubarn );
  fprintf( stdout, "------------------------------------\n" );
#endif 

  std::cout << "FillEvent start    ( " << evt_buffer.size()  << " evt)" << std::endl;
  if(flag_event == 1) FillEvent(evt_buffer);
  std::cout << "FillEvent finished ( " << evt_buffer.size()  << " evt)" << std::endl;

  return evt_buffer;
}


std::vector<SKSNSimSNEventVector> SKSNSimVectorSNGenerator::MakeEvent(double time, double nu_energy, int nReact, int nuType, double rate){
  // SKSNSimTools::DumpDebugMessage(Form(" MakeEvent time %.2g nuEne %.2g nReact %d nuType %d rate %.2g", time , nu_energy, nReact, nuType, rate));
  std::vector<SKSNSimSNEventVector> buffer;

  //double totcrsIBD[nuEneNBins] = {0.};
  double dRandTotEvts = randomgenerator->Poisson(rate);
  //if(time<0.005)std::cout << time << " " << nu_energy << " " << nReact << " " << nuType << " " << rate << std::endl; //nakanisi


  SKSNSimBinnedFluxModel &flux = dynamic_cast<SKSNSimBinnedFluxModel&>(*fluxmodels[0]); // TODO selectable flux
  const double nuEneBinSize = flux.GetBinWidthEne(0);
  const double tBinSize = flux.GetBinWidthTime(0);

  auto getRandomReal = [](double s, double e , TRandom& rng){ return rng.Uniform(s,e);};

  auto determinePosition = std::bind([](TRandom &rng)
  {
    const double rPositionRange = RINTK;
    const double hPositionRange = ZPINTK;
    //random inside the full tank (32.5kton)
    const double r2 = rng.Uniform(1.) * rPositionRange*rPositionRange ;
    const double r = std::sqrt( r2 );
    const double phi = rng.Uniform( 2. * M_PI);
    const double x = r * std::cos( phi );
    const double y = r * std::sin( phi );
    const double z = -hPositionRange + rng.Uniform( 2.*hPositionRange);

    return UtilVector3<double>(x,y,z);
  }, *randomgenerator);

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
      auto xyz = determinePosition();

      //SNEvtInfo evtInfo;
      //evtInfo.iEvt = dRandTotEvts;
      SKSNSimSNEventVector evtInfo;

      double rvtx [3];
      rvtx[0] = xyz.x;
      rvtx[1] = xyz.y;
      rvtx[2] = xyz.z;
      evtInfo.SetSNEvtInfo(nReact, tReact, nuType, nuEne, sn_dir, rvtx);

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
	int totGenNcNup=0, totGenNcNun=0, totGenNcNubarp=0, totGenNcNubarn=0;
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
    determineKinematics( xsecmodels, *randomgenerator, p, sn_dir);

    const auto rType = p.GetSNEvtInfoRType();
    int Reaction  = rType/10e4 - 1;
    int State_pre = rType/10e3;
    int Ex_state_pre = rType/10;
    int State = ((rType - (Reaction+1)*10e4)/10e3) - 1;
    int Ex_state = ((rType - State_pre*10e3)/10) - 1;
    int channel = (rType - Ex_state_pre*10) - 1;
    if(iSkip == 0) {

      if(rType == 0) totGenNuebarp++;
      if(rType == 1) totGenNueElastic++;
      if(rType == 2) totGenNuebarElastic++;
      if(rType == 3) totGenNuxElastic++;
      if(rType == 4) totGenNuxbarElastic++;
      if(rType>4){
        //if(Ex_state==30)std::cout << p.rType << " " << "nReact" << " " << nReact << " " << "Reaction" << " " << Reaction << " " << "State_pre" << " " << State_pre << " " << "State" << " " << State << " " << "Ex_state_pre" << " " << Ex_state_pre << " " << "Ex_state" << " " << Ex_state << " " << "channel" << " " << channel << std::endl; //nakanisi
        if(Reaction==0 && Ex_state!=29) totGenNueO++;
        if(Reaction==1 && Ex_state!=29) totGenNuebarO++;
        if(Reaction==0 && Ex_state==29) totGenNueOsub++;
        if(Reaction==1 && Ex_state==29) totGenNuebarOsub++; 
      }
    }
  }

  int totalNumOfGenEvts = ( totGenNuebarp + totGenNueElastic
      + totGenNuebarElastic + totGenNuxElastic + totGenNuxbarElastic
      + totGenNueO + totGenNuebarO 
      + totGenNueOsub + totGenNuebarOsub
      + totGenNcNup + totGenNcNun + totGenNcNubarp + totGenNcNubarn
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
  //fprintf( stdout, "   nue + o sub = %d\n", totGenNueOsub );
  //fprintf( stdout, "   nuebar + o sub = %d\n", totGenNuebarOsub );
  fprintf( stdout, "------------------------------------\n" );

}

void SKSNSimVectorSNGenerator::determineKinematics( std::map<XSECTYPE, std::shared_ptr<SKSNSimCrosssectionModel>> xsecmodels, TRandom &rng, SKSNSimSNEventVector &ev, const double snDir[])
{
  auto SQ = [](double x){return x*x;};
  const double nuEne = ev.GetSNEvtInfoNuEne();
  const auto pvect = nuEne * ev.GetSNEvtInfoNuDir();

  //number of particle emitted on deexcitation with CC reaction
  constexpr int numNtNueO[7] = {0, 1, 0, 2, 0, 0, 0};
  constexpr int numPtNueO[7] = {1, 1, 2, 1, 0, 0, 1};
  constexpr int numNtNuebarO[7] = {0, 1, 0, 2, 1, 0, 0};
  constexpr int numPtNuebarO[7] = {0, 0, 1, 0, 1, 1, 2};
  constexpr int numGmNuebarO[7] = {1, 0, 0, 0, 0, 0, 0};

  const UtilVector3<double> snDir_vec(snDir);
  int iSkip = 0;

  SKSNSimXSecIBDSV       &xsecibd         = dynamic_cast<SKSNSimXSecIBDSV&>(      *xsecmodels[mXSECIBD]);
  SKSNSimXSecNuElastic   &xsecnuela       = dynamic_cast<SKSNSimXSecNuElastic&>(  *xsecmodels[mXSECELASTIC]);
  SKSNSimXSecNuOxygen    &xsecnuoxygen    = dynamic_cast<SKSNSimXSecNuOxygen&>(   *xsecmodels[mXSECOXYGEN]);
  SKSNSimXSecNuOxygenSub &xsecnuoxygensub = dynamic_cast<SKSNSimXSecNuOxygenSub&>(*xsecmodels[mXSECOXYGENSUB]);

	double sn_theta = acos( snDir[2] );
	double sn_phi = atan2( snDir[1],  snDir[0] );
  const UtilMatrix3<double> Rmat( cos(sn_theta)*cos(sn_phi), -sin(sn_phi), sin(sn_theta)*cos(sn_phi),
                                  cos(sn_theta)*sin(sn_phi),  cos(sn_phi), sin(sn_theta)*sin(sn_phi),
                                             -sin(sn_theta),           0.,             cos(sn_theta));

  auto nReact = ev.GetSNEvtInfoRType();
  if( nReact == 0 ){ // nuebar + p -> e+ + n

    // Original neutrino
    const auto nuMomentum = pvect;
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
    determineAngleNuebarP( rng, xsecibd, nuEne, eEne, eTheta, ePhi );
    double amom = sqrt(SQ( eEne ) - SQ( Me ));
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

  } else {
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
    auto determineNeutMomentum = std::bind([](TRandom &rng)
        {
        //random reaction of neutron
        double phi = rng.Uniform(0., 2.*M_PI); 
        double theta = rng.Uniform(0., M_PI); // TODO this should be dCosTheta? This is based on original code
        return UtilVector3<double>(theta, phi);
        }, rng);

    if(numNtNueO[channel]!=0 || numNtNuebarO[channel]!=0 || numGmNuebarO[channel]!=0){
      if(Reaction==0){
        int i_nucre = 0;
        if(Ex_state!=29){
          for(int i=0;i<numNtNueO[channel];i++){
            // Neutron
            i_nucre++;
            const auto neutronMom = sqrt(SQ(0.5+Mn) - SQ(Mn)) * determineNeutMomentum();
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
              auto neutronMom = sqrt(SQ(0.5+Mn)-SQ(Mn)) * determineNeutMomentum();
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
            const auto gammaMom = gammaEne * determineNeutMomentum();
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

