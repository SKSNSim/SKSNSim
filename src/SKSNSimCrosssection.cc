/**************************************
 * File: SKSNSimCrosssection.cc
 * Description:
 * Cross section interface
 *************************************/

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <pdg_codes.h>
#include "SKSNSimCrosssection.hh"
#include "SKSNSimConstant.hh"

using namespace SKSNSimPhysConst;

namespace SKSNSimCrosssection {
  double CalcIBDEnuFromEpos ( const double Epos /* MeV */, const double cosTheta ){
    /* Just copied from VectGenNuCrosssection::getEnu */ // TODO check if it is correct
    constexpr double delta = Mn-Mp;
    constexpr double d = (delta*delta - Me*Me)/(2*Mp);
    double pe=sqrt(Epos*Epos-Me*Me);

    return((Epos+delta+d)/(1.-(Epos-pe*cosTheta)/Mp));
  }
}

std::pair<double,double> SKSNSimXSecIBDSV::GetDiffCrosssection(double enu /* MeV */ , double costheta) const {
  /*********************************************************************/
  /// Copied from VectGenNuCrosssection::DcsNuebP_SV
  /// (purpose)
  ///     Calculate IBD cross section
  /// (input)
  ///     ENU, COSTHETA
  /// (output)
  ///     std::pair<SIGM,EP>
  /// 2019.04 production by H.Ito
  /// 2019.05 inport to slmkb8ibd.F by H.Ito
  /// 2020.05 rewrite to C++ by H.MHarada
  /*********************************************************************/


  // Computations using masses   
  constexpr double MnMn = Mn*Mn;
  constexpr double MeMe = Me*Me;
  constexpr double MpMp = Mp*Mp;
  constexpr double ave_nucleon_mass = (Mn + Mp)/2.;
  constexpr double Delta_CM = (MnMn - MpMp - MeMe)/(2.*Mp);

  //      Fermi coupling constant
  constexpr double HBAR_C2 = 0.3893793656e-21;// (MeV^2 * cm^2) ... ( plank-constant * photon-velocity )^2
  constexpr double GF_GeV = 1.1663787e-11;// (MeV^-2)
  constexpr double GF = GF_GeV*GF_GeV * HBAR_C2;// *1e20*1e20;//(MeV^-2 * cm^2)

  constexpr double ANO_NU = (2.7928473446-1.)+1.9130427;// sentence btw formula (7) and (8)
  constexpr double ALPHA = 1. / 137.035999139;

  //MESH_ES = 100;
  // constexpr double RHO_NC = 1.0126;


  //constexpr double IBD_THR  = ((N_MASS+E_MASS)*(N_MASS+E_MASS)- P_MASS*P_MASS)/ (2.*P_MASS);
  //constexpr double ETHRE = 0.;

  // Positron energy:                                                                                                                                                                    
  const double epsilon = enu/Mp;// formula (8)
  const double kappa = std::pow(1.+epsilon,2) - std::pow(epsilon*costheta,2);// in sentence below formula (21)
  const double Epo = ((enu-Delta_CM)*(1.+epsilon)+epsilon*costheta*sqrt(pow(enu-Delta_CM,2)-MeMe*kappa))/kappa;// formula (21)

  // Parameters
  const double Ppo = sqrt(Epo*Epo-MeMe);// formula (21)
  const double dE_dCosT = Ppo*epsilon/(1.+epsilon*(1.-costheta*Epo/Ppo));// trans. formula (20)

  // Belows are written in formula (3) - (11)
  const double s_minus_u = 2.*Mp*(enu+Epo)-MeMe;// above of formula (11)
  const double t = MnMn-MpMp-2.*Mp*(enu-Epo);// above of formula (11)

  // formula (7)      
  const double x  = t /(4.*pow(ave_nucleon_mass,2.)); // part of numerator of 1st formula (7) 
  const double y  = 1. - t/710000.;// right side of denominator of 1st formula (7) 
  const double z  = 1. - t/1030000.;// denominator of 2nd formula (7)
  //y  = 1. - t/710**2;// right side of denominator of 1st formula (7) 
  //z  = 1. - t/1030**2;// denominator of 2nd formula (7)
  const double f1 = (1. - (1.+ANO_NU)*x)/((1.-x)*(y*y));// 1st formula (7)
  const double f2 = ANO_NU/((1.-x)*(y*y));// 1st formula (7)
  const double g1 = -1.27/z*z;// 2nd formula (7)
  const double g2 = 2.*g1*pow(ave_nucleon_mass,2.)/(Mpi*Mpi-t);

  //bolow A,B,C calculation is formula (6)      
  const double A = 1./16.*((t-MeMe)
            *(4.*pow(f1,2.)*(4.*pow(ave_nucleon_mass,2.)+t+MeMe)
                + 4.*pow(g1,2.)*(-4*pow(ave_nucleon_mass,2.)+t+MeMe)
                + pow(f2,2.)*(t*t/pow(ave_nucleon_mass,2.)+4.*t+4.*MeMe)
                + 4*MeMe*t*pow(g2,2.)/pow(ave_nucleon_mass,2.)
                + 8.*f1*f2*(2*t+MeMe)
                + 16*MeMe*g1*g2)
            - pow(DeltaM,2.)
            *((4.*pow(f1,2.)+t*pow(f2,2.)/pow(ave_nucleon_mass,2.))*(4.*pow(ave_nucleon_mass,2.)+t-MeMe)
                + 4.*pow(g1,2.)*(4.*pow(ave_nucleon_mass,2.)-t+MeMe)
                + 4.*MeMe*pow(g2,2.)*(t-MeMe)/pow(ave_nucleon_mass,2.)
                + 8.*f1*f2*(2.*t-MeMe)
                + 16.*MeMe*g1*g2)
            - 32.*MeMe*ave_nucleon_mass*DeltaM*g1*(f1+f2));

  const double B = 1./16.*(16.*t*g1*(f1+f2) + 4.*MeMe*DeltaM*(pow(f2,2.)+f1*f2+2.*g1*g2)/ave_nucleon_mass);

  const double C = 1./16.*(4.*(pow(f1,2)+pow(g1,2)) - t*pow(f2,2.)/pow(ave_nucleon_mass,2.));

  const double abs_M_squared = A - B*s_minus_u + C*pow(s_minus_u,2.);// formula (5)
  const double rad_correction = ALPHA/M_PI*(6.00352 + 3./2.*log(Mp/(2.*Epo)) + 1.2*pow(Me/Epo,1.5));// formula (14)        

  const double dsigma_dt = GF*pow(costheta_cabibo,2.)/(2*M_PI*std::pow(2.*Mp*enu,2.))*abs_M_squared;// formula (3) using s-Mp^2 = 2Mp*Enu(above of formula (11))
  const double dsigma_dE = 2.*Mp*dsigma_dt;//formula (11)
  double dcs = dE_dCosT*dsigma_dE*(1.+rad_correction);// result 
  //IBD_SIG0 = GF * COSTHETA_C**2 /FUNC_PI * (2. * P_MASS) / (8. * P_MASS**2);
  //SIGM = dE_dCosT * IBD_SIG0 / ENU**2 * abs_M_squared * (1. + rad_correction)  

  if ( dcs < 0 ) dcs = 0.;
#ifdef DEBUG
  std::cout << "In SKSNSimXSecIBDSV: (sigm, Epo) = ( " << dcs << ", " << Epo << " )" << std::endl;
#endif
  return std::make_pair(dcs, Epo);
}

double SKSNSimXSecIBDSV::GetCrosssection(double enu) const
{
  /*
     Total cross section of nu_e_bar + p --> e^+  + n  interaction
     using the precise calculation by Strumia and Vissani
     enu : Neutrino energy in MeV
     (Epo : Positoron total energy in MeV, should be later)
     return : Total cross section (cm^2)
  */

  double totcsnuebp_SV = 0;
  const double Ee = enu - DeltaM;
  if(Ee <= Me) return 0.0;

  constexpr double COSTHETARANGELOW  = -1.;
  constexpr double COSTHETARANGEHIGH =  1.;
  constexpr int    NBIN = 200;
  constexpr double BINWIDTH = (COSTHETARANGEHIGH - COSTHETARANGELOW)/double(NBIN);
  constexpr double HALFBINWIDTH = BINWIDTH / 2.0;
  constexpr double COSTHETARANGELOWINCLBIAS = COSTHETARANGELOW + HALFBINWIDTH;
  for(int b=0; b<NBIN; b++){
    double costheta = COSTHETARANGELOWINCLBIAS + BINWIDTH * double(b);
    double dcs0, Epo0;
    auto result = GetDiffCrosssection(enu, costheta);
    totcsnuebp_SV += result.first;
  }
  totcsnuebp_SV *= BINWIDTH;

  //cout << totcsnuebp_SV << endl;
  return totcsnuebp_SV;
}

std::pair<double,double> SKSNSimXSecIBDVB::GetDiffCrosssection(double enu, double costheta)
 const {
  /*
    get differential cross section of nue_bar P and positron energy
         enu : Neutrino energy in MeV
    costheta : 
         Epo : Positron total energy in MeV
         dcs : Differenctial cross section d sigma / d cos(theta) (cm^2)

    "The angular distribution of the reaction nu_e_bar + p -> e+ + n"
     P.Vogel and J.F.Beacom [arXiv:hep-ph/9903554 1 Apr 1999]
  */

  const double Epo0 = enu - DeltaM;
  constexpr double MeMe = Me*Me;
  constexpr double MpMp = Mp*Mp;
  constexpr double Delta_R_inner = 0.024; // inner radiative correction
  constexpr double f     = 1.;
  constexpr double g     = 1.26; // axial-vector coupling constant

  if(Epo0 < Me) return std::make_pair(0.0, 0.0);

  const double Pe0 = sqrt(pow(Epo0, 2.) - MeMe); // positron momemtum

  // formula(9)
  const double sigma0 = pow(Gf,2.) * pow(costheta_cabibo,2.) / PI * (1. + Delta_R_inner) / pow(HBARC,2.) * 1.0E-26;

  // formula(13)
  constexpr double ave_nucleon_mass = (Mp + Mn)/2.; // average nucleon mass
  const double v_e_0    = Pe0 / Epo0; // velocity
  const double y        = sqrt((pow(DeltaM,2.) - MeMe)/2.);

  const double Epo = Epo0 * (1.- enu / ave_nucleon_mass * (1.- v_e_0 *costheta)) - pow(y,2.) / ave_nucleon_mass;

  if(Epo < Me) return std::make_pair(0.0, 0.0);


  // formula(15)
  constexpr double f2   = 3.706; // nu_p - nu_n (magnetic moment)
  const double gamma = 2*(f + f2) * g * ((2*Epo0 + DeltaM) * (1. - v_e_0 * costheta) - MeMe / Epo0) + (pow(f,2.) + pow(g,2.)) * (DeltaM * (1. + v_e_0 * costheta) + MeMe / Epo0) + (pow(f,2.) + 3 * pow(g,2.)) * ((Epo0 + DeltaM)*(1. - 1./v_e_0 * costheta) - DeltaM) + (pow(f,2.) - pow(g,2.)) * ((Epo0 + DeltaM) * (1. - 1./v_e_0 * costheta) - DeltaM) * v_e_0 * costheta;

  // formula(14)
  const double Pe1 = sqrt(pow(Epo, 2.) - MeMe); // positron momemtum
  const double v_e_1 = Pe1 / Epo; // velocity

  const double dcs = sigma0 / 2. * ((pow(f,2.) + 3 * pow(g,2.)) + (pow(f,2.) - pow(g,2.)) * v_e_1 * costheta) * Epo * Pe1 - sigma0 / 2 *(gamma / ave_nucleon_mass) * Epo0 * Pe0;
  //cout << dcs << " " << sigma0 << " " << gamma << " " << v_e_1 << " " << costheta << endl;

  return std::make_pair(dcs, Epo);
}

double SKSNSimXSecIBDVB::GetCrosssection(double enu) const
{
  /*
    Total cross section of nu_e_bar + p --> e^+  + n  interaction
    using the precise calculation by Vogel and Beacom
        enu : Neutrino energy in MeV
	(Epo : Positoron total energy in MeV, should be later)
     return : Total cross section (cm^2)
  */

  double totcsnuebp_vb = 0;
  const double Ee = enu - DeltaM;
  if(Ee < 3.0 /* MeV */) return 0.0;

  constexpr double COSTHETARANGELOW  = -1.;
  constexpr double COSTHETARANGEHIGH =  1.;
  constexpr int    NBIN = 200;
  constexpr double BINWIDTH = (COSTHETARANGEHIGH - COSTHETARANGELOW)/double(NBIN);
  constexpr double HALFBINWIDTH = BINWIDTH / 2.0;
  constexpr double COSTHETARANGELOWINCLBIAS = COSTHETARANGELOW + HALFBINWIDTH;
  for(int b=0; b<NBIN; b++){
    double costheta = COSTHETARANGELOWINCLBIAS + BINWIDTH * double(b);
    double dcs0, Epo0;
    auto result = GetDiffCrosssection(enu, costheta);
    totcsnuebp_vb += result.first;
  }
  totcsnuebp_vb *= BINWIDTH;

  //cout << totcsnuebp_SV << endl;
  return totcsnuebp_vb;
}

double SKSNSimXSecNuElastic::GetCrosssection(double enu, int ipart, FLAGETHR flag) const
{
  /*
    Total cross section of nu + e --> nu + e  interaction
    using the precise calculation by Bahcall refer from sollib/sl_elctot_rad.F
        enu : Neutrino energy in MeV
     return : Total cross section (cm^2)
  */

  if(!((ipart == PDG_ELECTRON_NEUTRINO) || (ipart == - PDG_ELECTRON_NEUTRINO) || (ipart == PDG_MUON_NEUTRINO) || (ipart == -PDG_MUON_NEUTRINO))){
    std::cerr << "Not support particle code " << ipart << std::endl;
    exit(1);
  }

  double x = 0.;
  if(enu <= nuElaEneMin || enu >= nuElaEneMax) return x;

  if(flag == ETHRON) { // should apply Eth
    constexpr double t_min = eEneThr - Me;
    const double t_max = 2.*enu*enu/(Me+2.*enu);

    if(t_min > t_max){
      x = 0.;
      return x;
    }

    constexpr int istep = 1000;
    const double dstep = (t_max - t_min) / double(istep);
    x=0.;
    for (int i = 0; i<istep; i++){
      double E = (t_min + dstep/2.) + dstep * double(i) + Me;
      switch (ipart)
      {
        case  PDG_ELECTRON_NEUTRINO: x += sl_nue_dif_rad_(&enu, &E) * dstep; break;
        case -PDG_ELECTRON_NEUTRINO: x += sl_neb_dif_rad_(&enu, &E) * dstep; break;
        case  PDG_MUON_NEUTRINO:     x += sl_num_dif_rad_(&enu, &E) * dstep; break;
        case -PDG_MUON_NEUTRINO:     x += sl_nmb_dif_rad_(&enu, &E) * dstep; break;
        default: std::cerr << "Should not appear this message: " << __FILE__ << " L:" << __LINE__ << std::endl; break;
      }
    }
  }
  else if(flag == ETHROFF) { // should not apply Eth
    double e1, e2, c1=0., c2=0.;

    const int iene = (int)(enu / nuElaEneBinSize);
    fCsElaTree->GetEntry(iene-1);
    e1 = NuElaEnergy;
    if(ipart ==  12) c1 = CsElaNue;
    if(ipart == -12) c1 = CsElaNeb;
    if(ipart ==  14) c1 = CsElaNux;
    if(ipart == -14) c1 = CsElaNxb;

    fCsElaTree->GetEntry(iene);
    e2 = NuElaEnergy;
    if(ipart ==  12) c2 = CsElaNue;
    if(ipart == -12) c2 = CsElaNeb;
    if(ipart ==  14) c2 = CsElaNux;
    if(ipart == -14) c2 = CsElaNxb;

    x = (c2-c1) / (e2-e1) * (enu - e1) + c1;
    //std::cout << enu << " " << iene << " " << e1 << " " << e2 << " " << c1 << " " << c2 << " " << x << std::endl;
  }

  return x;

}

void SKSNSimXSecNuElastic::OpenCsElaFile(){
  fCsElaFile.reset(new TFile(CSELAFILENAME.c_str(), "READ"));
  fCsElaTree = (TTree*)fCsElaFile->Get("nuela");
  fCsElaTree->SetBranchAddress("nuene", &NuElaEnergy);
  fCsElaTree->SetBranchAddress("cnue", &CsElaNue);
  fCsElaTree->SetBranchAddress("cnux", &CsElaNux);
  fCsElaTree->SetBranchAddress("cneb", &CsElaNeb);
  fCsElaTree->SetBranchAddress("cnxb", &CsElaNxb);
}

void SKSNSimXSecNuElastic::CloseCsElaFile(){
  delete fCsElaTree;
  if( fCsElaFile != NULL)
    fCsElaFile->Close();

}

std::pair<double,double> SKSNSimXSecNuElastic::GetDiffCrosssection(double e, double r) const{
  std::cerr << "Not supported: " << __FILE__ << " GetDiffCrosssection(...) for NuElastic XSec model" << std::endl;
  return std::make_pair(0.0,0.0);
}

double SKSNSimXSecNuElastic::CalcElectronTotEnergy(const double nuEne, const double cost)
{
  auto SQ = [] (const double &x){return x*x;};
	double alpha = Me / nuEne;
	double e; 
	if( cost > 0. ){
		e = 2. * alpha * SQ( cost ) 
			/ ( SQ( 1. + alpha ) - SQ( cost ) )
			* nuEne 
			+ Me; //total energy
	}else{
		e = 0.;
	}
	return e;
}

double SKSNSimXSecNuElastic::CalcDeEneDCost(const double nuEne, const double cost)
{
  auto SQ = [] (const double &x){return x*x;};
	double alpha = Me / nuEne;
	double p = nuEne * 4. * alpha * cost * SQ( 1 + alpha )
		/ SQ( SQ( 1. + alpha ) - SQ( cost ) );
	return p;
}
double SKSNSimXSecNuElastic::CalcCosThr(const double nuEne, const double eEth){
	double alpha = Me / nuEne;
	double y = ( eEth - Me ) / nuEne;

	double cosTth = sqrt( y / ( y + 2. * alpha ) ) * ( 1. + alpha );

	return cosTth;
}

void SKSNSimXSecNuOxygen::LoadFile(){
  for(int i = 0; i < NTYPE; i++){
    for(int j = 0; j < NIXSTATE; j++){
      LoadFile( convToINISTATE(i,j));
    }
  }
}
void SKSNSimXSecNuOxygen::LoadFile(INISTATE ini){
  const static std::string rctn[NTYPE] = {"nue", "neb"};
  if (!isOpened[ini]) {
    std::cout << "new file open (num, ix, isOpen): " << std::get<0>(ini) << " " << std::get<1>(ini) << " " << isOpened[ini] << std::endl;
    const std::string target = Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[std::get<0>(ini)].c_str(),std::get<1>(ini));
    std::ifstream ifs(target);

    if(!ifs){
      std::cerr << "file does not exist" << std::endl;
      std::cerr << "file name is" << " " << target << std::endl;
      exit(EXIT_FAILURE);
    }

    double tmp, tmp_e0;
    double tmp_num, tmp_ex, tmp_rec, tmp_pro, tmp_sum;
    double tmp_cs[NCHANNEL];
    int index = 0;
    while(ifs>>tmp>>tmp_e0){
      nuene[ini].push_back(tmp_e0);
      for(int j=0;j<GetNumEx(std::get<1>(ini));j++){
        ifs >>tmp_num>>tmp_ex>>tmp_rec>>tmp_pro>>tmp_sum>>tmp_cs[0]>>tmp_cs[1]>>tmp_cs[2]>>tmp_cs[3]>>tmp_cs[4]>>tmp_cs[5]>>tmp_cs[6];
        index++;
        exEne[ini].push_back(tmp_ex);
        rec[ini].push_back(tmp_rec);
        pro_energy[ini].push_back(tmp_pro);
        for(int k = 0; k < NCHANNEL; k++){
          INIFINSTATE inifin = convToINIFINSTATE(ini, j, k);
          crs[inifin].push_back(tmp_cs[k]);
        }
      }
    }
    fileSize[ini] = index;
    //std::cout << fileSize[num][ix] << std::endl;
    //std::cout << num << " " << ix << " " << index << std::endl;
    //std::cout << exEne[num][ix].at(0) << " " << exEne[num][ix].at(1) << " " << exEne[num][ix].at(2) << std::endl;
    //
    isOpened[ini] = true;
    ifs.close();
  }
}

void SKSNSimXSecNuOxygen::InitializeTable(){
  isOpened.clear();
  rec.clear();
  pro_energy.clear();
  fileSize.clear();
  exNum.clear();
  exEne.clear();
  nuene.clear();
  crs.clear();
  for(int t = 0; t < NTYPE; t++){
    exNum[convToINISTATE(t,0)] = 3;
    exNum[convToINISTATE(t,1)] = 15;
    exNum[convToINISTATE(t,2)] = 8;
    exNum[convToINISTATE(t,3)] = 1;
    exNum[convToINISTATE(t,4)] = 6;
    for(int ix = 0; ix < NIXSTATE; ix++){
      INISTATE ini = convToINISTATE(t,ix);
      isOpened[ini] = false;
      rec[ini].clear();
      pro_energy[ini].clear();
      exEne[ini].clear();
      nuene.clear();
      fileSize[ini] = 0;
      for(int ex = 0; ex < NEXSTATE; ex++){
        for(int ch = 0; ch < NCHANNEL; ch++){
          INIFINSTATE inifin = convToINIFINSTATE(ini,ex,ch);
          crs[inifin].clear();
        }
      }
    }
  }
}

double SKSNSimXSecNuOxygen::GetCrosssection(double enu, INIFINSTATE inifin) const
{
  /*
   * Total cross section of nu_e + O --> e^- + X, nu_e_bar + O --> e^+ + X interaction
   */

  double totcso = 0;

  const INISTATE ini = GetINISTATE(inifin);
  const int ex = std::get<2>(inifin);
  const int ch = std::get<3>(inifin);
  double rec_energy = 0.0;
  auto findBin = [ini] (double e, const std::vector<double> &enebin){
    for(int i = 0; i < 20; i++){
      if( i > 0 && e < enebin.at(i) && e >= enebin.at(i-1)) return i;
    }
    return -1;
  };

  int ienebin = findBin(enu, nuene.at(ini)); 
  if( ienebin == -1 ) return .0;
  rec_energy = enu - exEne.at(ini).at(ex);
  if(rec_energy > Me && rec_energy>eEneThr){
    totcso = (((crs.at(inifin).at(ienebin)-crs.at(inifin).at(ienebin-1))/(nuene.at(ini).at(ienebin)-nuene.at(ini).at(ienebin-1)))*enu + crs.at(inifin).at(ienebin-1) - ((crs.at(inifin).at(ienebin)-crs.at(inifin).at(ienebin-1))/(nuene.at(ini).at(ienebin)-nuene.at(ini).at(ienebin-1)))*nuene.at(ini).at(ienebin-1)) * 1.0e-26;
  }
  else
    totcso = 0.;
  return totcso;
}

std::pair<double,double> SKSNSimXSecNuOxygen::GetDiffCrosssection(double e, double r) const
{
  std::cerr << "Not supported: " << __FILE__ << " GetDiffCrosssection(...) for NuOxygen XSec model" << std::endl;
  return std::make_pair(0.0,0.0);
}

double SKSNSimXSecNuOxygen::OxigFuncAngleRecCC(int num, int ix, int ex, int ch, double enu, double cos)
{
  /*
   * calculate probability of the direction the electron is emitted after the nu_e + O charged current reaction.
   */

  double costheta = 0;

  std::ifstream ifs[2][5];
  std::string rctn[2] = {"nue", "neb"};

  ifs[num][ix].open(Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(),ix));

  if(!ifs[num][ix]){
    std::cerr << "file does not exist" << std::endl;
    std::cerr << "file name is" << " " << Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(),ix) << std::endl;
    exit(EXIT_FAILURE);
  }

  int fileSize[2][5] = {{0}};
  int exNum[2][5] = {{3,15,8,1,6}};
  int num_ex = 0;
  double nuEne = 0.;
  std::vector<double> Enu[2][5];
  std::vector<double> totcrs[2][5];
  std::vector<double> Ee[2][5];
  int index = 0;
  double tmp_ene, tmp_crs;
  std::vector<double> nuene[2][5];
  double tmp_e0;
  double tmp;
  double tmp_cs0, tmp_cs1, tmp_cs2, tmp_cs3, tmp_cs4, tmp_cs5, tmp_cs6, tmp_sum, tmp_rec, tmp_pro, tmp_ex, tmp_num;
  std::vector<double> exEne[2][5];
  std::vector<double> rec[2][5];
  std::vector<double> pro_energy[2][5];
  if(ix==0) num_ex = 3;
  else if(ix==1) num_ex = 15;
  else if(ix==2) num_ex = 8;
  else if(ix==3) num_ex = 1;
  else if(ix==4) num_ex = 16;
  //std::cout << ix << " " << "num_ex is" << " " << num_ex << std::endl;
  while(ifs[num][ix]>>tmp>>tmp_e0){
    nuene[num][ix].push_back(tmp_e0);
    for(int j=0;j<num_ex;j++){
      ifs[num][ix]>>tmp_num>>tmp_ex>>tmp_rec>>tmp_pro>>tmp_sum>>tmp_cs0>>tmp_cs1>>tmp_cs2>>tmp_cs3>>tmp_cs4>>tmp_cs5>>tmp_cs6;
      index++;
      exEne[num][ix].push_back(tmp_ex);
      rec[num][ix].push_back(tmp_rec);
      pro_energy[num][ix].push_back(tmp_pro);
    }
  }
  fileSize[num][ix] = index;

  constexpr double me = Me; // MeV

  if(ch!=8){
    double rec_energy[16] = {0.};
    rec_energy[ex] = enu - exEne[num][ix].at(ex); //energy [MeV] of recoil electron or positron

    double a = 1. + (rec_energy[ex]/25.)*(rec_energy[ex]/25.)*(rec_energy[ex]/25.)*(rec_energy[ex]/25.);
    double b = 3. + (rec_energy[ex]/25.)*(rec_energy[ex]/25.)*(rec_energy[ex]/25.)*(rec_energy[ex]/25.);

    costheta = 0.5 * (1.-(a/b)*cos);

    //return costheta;
  }

  if(ch==8){
    double recEnergy = 0.;
    if(num==0){
      recEnergy = enu - 15.4;
    }
    if(num==1){
      recEnergy = enu - 11.4;
    }
    double a = 1. + (recEnergy/25.)*(recEnergy/25.)*(recEnergy/25.)*(recEnergy/25.);
    double b = 3. + (recEnergy/25.)*(recEnergy/25.)*(recEnergy/25.)*(recEnergy/25.);

    costheta = 0.5 * (1.-(a/b)*cos);

  }
  return costheta;
}

double SKSNSimXSecNuOxygen::OxigFuncRecEneCC(int num, int ix, int ex, int ch, double enu)
{
  /*
   * calculate recoil e- or e+ energy 
   */

  double recEnergy = 0;

  if(ch!=8){
    std::ifstream ifs[2][5];
    std::string rctn[2] = {"nue", "neb"};

    ifs[num][ix].open(Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(),ix));

    if(!ifs[num][ix]){
      std::cerr << "file does not exist" << std::endl;
      std::cerr << "file name is" << " " << Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(), ix) << std::endl;
      exit(EXIT_FAILURE);
    }

    int fileSize[2][5] = {{0}};
    int exNum[2][5] = {{3,15,8,1,6}};
    int num_ex = 0;
    double nuEne = 0.;
    std::vector<double> Enu[2][5];
    std::vector<double> Ee[2][5];
    int index = 0;
    double tmp_ene, tmp_crs;
    std::vector<double> nuene[2][5];
    double tmp_e0;
    double tmp;
    double tmp_cs0, tmp_cs1, tmp_cs2, tmp_cs3, tmp_cs4, tmp_cs5, tmp_cs6, tmp_sum, tmp_rec, tmp_pro, tmp_ex, tmp_num;
    std::vector<double> exEne[2][5];
    std::vector<double> rec[2][5];
    std::vector<double> pro_energy[2][5];
    if     (ix==0) num_ex = 3;
    else if(ix==1) num_ex = 15;
    else if(ix==2) num_ex = 8;
    else if(ix==3) num_ex = 1;
    else if(ix==4) num_ex = 16;
    while(ifs[num][ix]>>tmp>>tmp_e0){
      nuene[num][ix].push_back(tmp_e0);
      for(int j=0;j<num_ex;j++){
        ifs[num][ix]>>tmp_num>>tmp_ex>>tmp_rec>>tmp_pro>>tmp_sum>>tmp_cs0>>tmp_cs1>>tmp_cs2>>tmp_cs3>>tmp_cs4>>tmp_cs5>>tmp_cs6;
        index++;
        exEne[num][ix].push_back(tmp_ex);
        rec[num][ix].push_back(tmp_rec);
        pro_energy[num][ix].push_back(tmp_pro);
      }
    }
    fileSize[num][ix] = index;

    constexpr double me = Me; // MeV

    double rec_energy[16] = {0.};
    rec_energy[ex] = enu - exEne[num][ix].at(ex); //energy [MeV] of recoil electron or positron

    recEnergy = rec_energy[ex];
  }

  if(ch==8){
    if(num==0){
      if(enu > 15.4){
        recEnergy = enu - 15.4;
      }
      else recEnergy = 0.;
    }
    if(num==1){
      if(enu > 11.4){
        recEnergy = enu - 11.4;
      }
      else recEnergy = 0.;
    }
  }

  return recEnergy;
}

void SKSNSimXSecNuOxygenNC::LoadFile(INISTATE ini){
  const static std::string rctn[NTYPE] = {"N", "O"};
  if (!isOpened[ini]) {
    std::cout << "new file open (num, isOpen): " << std::get<0>(ini) << " " << isOpened[ini] << std::endl;
    const std::string target = Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crossNC_ex%s.dat",rctn[std::get<0>(ini)].c_str());
    std::ifstream ifs(target);

    if(!ifs){
      std::cerr << "file does not exist" << std::endl;
      std::cerr << "file name is" << " " << target << std::endl;
      exit(EXIT_FAILURE);
    }

    double tmp, tmp_e0;
    double tmp_ex;
    double tmp_cs;
    int index = 0;
    while(ifs>>tmp_e0){
      nuene[ini].push_back(tmp_e0);
      for(int j=0;j<GetNumEx(std::get<0>(ini));j++){
        INIFINSTATE inifin = convToINIFINSTATE(ini, j);
        ifs >>tmp_ex>>tmp_cs;
        index++;
        exEne[ini].push_back(tmp_ex);
        crs[inifin].push_back(tmp_cs);
      }
    }
    fileSize[ini] = index;

    isOpened[ini] = true;
    ifs.close();
  }
}

void SKSNSimXSecNuOxygenNC::LoadFile() {
  for( size_t t = 0; t < NTYPE; t++) LoadFile({t});
}

void SKSNSimXSecNuOxygenNC::InitializeTable(){
  isOpened.clear();
  fileSize.clear();
  exEne.clear();
  nuene.clear();
  crs.clear();
  for(int t = 0; t < NTYPE; t++){
    INISTATE ini = convToINISTATE(t);
    isOpened[ini] = false;
    fileSize[ini] = 0;
    for(int ch = 0; ch < GetNumEx(t); ch++){
      INIFINSTATE inifin = convToINIFINSTATE(ini,ch);
      crs[inifin].clear();
    }
  }
}

double SKSNSimXSecNuOxygenNC::GetCrosssection(double e, INIFINSTATE inifin) const {
  double totcso = 0.;

  const INISTATE ini = GetINISTATE(inifin);
  auto findBin = [] (double e, const std::vector<double> &enebin){
    for(size_t i = 1; i < enebin.size(); i++){
      if(e < enebin.at(i) && e >= enebin.at(i-1)) return (int)i;
    }
    return -1;
  };
  const int e_bin = findBin(e, nuene.at(ini));

  if( e_bin >= 0 && e>exEne.at(ini).at(e_bin)){
    totcso = ((((crs.at(inifin).at(e_bin)-crs.at(inifin).at(e_bin-1))/(nuene.at(ini).at(e_bin)-nuene.at(ini).at(e_bin-1)))*e + crs.at(inifin).at(e_bin-1) - ((crs.at(inifin).at(e_bin)-crs.at(inifin).at(e_bin-1))/(nuene.at(ini).at(e_bin)-nuene.at(ini).at(e_bin-1)))*nuene.at(ini).at(e_bin-1))) * 1.0e-42;
  }
  return totcso;
}

std::pair<double,double> SKSNSimXSecNuOxygenNC::GetDiffCrosssection(double e, double r) const
{
  std::cerr << "Not supported: " << __FILE__ << " GetDiffCrosssection(...) for NuOxygenNC XSec model" << std::endl;
  return std::make_pair(0.0,0.0);
}

void SKSNSimXSecNuOxygenSub::LoadFile(){
  for(int i = 0; i < NTYPE; i++){
    for(int j = 0; j < NIXSTATE; j++){
      LoadFile(convToINISTATE(i,j));
    }
  }
}
void SKSNSimXSecNuOxygenSub::LoadFile(INISTATE ini){
  const static std::string rctn[NTYPE] = {"nue", "neb"};
  if (!isOpened[ini]) {
    std::cout << "new file open (num, ix, isOpen): " << std::get<0>(ini) << " " << std::get<1>(ini) << " " << isOpened[ini] << std::endl;
    const std::string target = Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/sub3-1/dat/sub%s%dstate.dat",rctn[std::get<0>(ini)].c_str(),std::get<1>(ini));
    std::ifstream ifs(target);

    if(!ifs){
      std::cerr << "file does not exist" << std::endl;
      std::cerr << "file name is" << " " << target << std::endl;
      exit(EXIT_FAILURE);
    }

    double tmp, tmp_e0;
    //double tmp_num, tmp_ex, tmp_rec, tmp_pro, tmp_sum;
    double tmp_cs[NCHANNEL];
    int index = 0;
    while(!ifs.eof()){
      ifs >>tmp_e0>>
        tmp_cs[0]>>tmp_cs[1]>>tmp_cs[2]>>tmp_cs[3]>>tmp_cs[4]
        >>tmp_cs[5]>>tmp_cs[6]>>tmp_cs[7]>>tmp_cs[8]>>tmp_cs[9]
        >>tmp_cs[10]>>tmp_cs[11]>>tmp_cs[12]>>tmp_cs[13]>>tmp_cs[14]
        >>tmp_cs[15]>>tmp_cs[16]>>tmp_cs[17]>>tmp_cs[18]>>tmp_cs[19]
        >>tmp_cs[20]>>tmp_cs[21]>>tmp_cs[22]>>tmp_cs[23]>>tmp_cs[24]
        >>tmp_cs[25]>>tmp_cs[26]>>tmp_cs[27]>>tmp_cs[28]>>tmp_cs[29]
        >>tmp_cs[30]>>tmp_cs[31];

      index++;
      Enu[ini].push_back(tmp_e0);
      for(int k = 0; k < NCHANNEL; k++){
        INIFINSTATE inifin = convToINIFINSTATE(ini, k);
        crs[inifin].push_back(tmp_cs[k]);
      }
    }
    fileSize[ini] = index-1;
    //std::cout << fileSize[num][ix] << std::endl;
    //std::cout << num << " " << ix << " " << index << std::endl;
    //std::cout << exEne[num][ix].at(0) << " " << exEne[num][ix].at(1) << " " << exEne[num][ix].at(2) << std::endl;
    //
    isOpened[ini] = true;
    ifs.close();
  }
}

void SKSNSimXSecNuOxygenSub::InitializeTable(){
  isOpened.clear();
  fileSize.clear();
  Enu.clear();
  crs.clear();
  for(int t = 0; t < NTYPE; t++){
    for(int ix = 0; ix < NIXSTATE; ix++){
      INISTATE ini = convToINISTATE(t,ix);
      isOpened[ini] = false;
      Enu[ini].clear();
      fileSize[ini] = 0;
      for(int ch = 0; ch < NCHANNEL; ch++){
        INIFINSTATE inifin = convToINIFINSTATE(ini,ch);
        crs[inifin].clear();
      }
    }
  }
}

double SKSNSimXSecNuOxygenSub::GetCrosssection(double enu, INIFINSTATE inifin) const
{
  /*
   * Total cross section of nu_e + O --> e^- + X, nu_e_bar + O --> e^+ + X interaction
   */

  double totcso = 0;

  const INISTATE ini = GetINISTATE(inifin);
  const int ch = std::get<2>(inifin);
  double rec_energy = 0.0;
  auto findBin = [ini] (double e, const std::vector<double> &enebin){
    for(size_t i = 0; i < enebin.size(); i++){
      if( i > 0 && e < enebin.at(i) && e >= enebin.at(i-1)) return (int)i;
    }
    return -1;
  };

  int ienebin = findBin(enu, Enu.at(ini)); 
  if( ienebin == -1 ) return .0;
  switch( std::get<0>(ini)){ // RCN
    case 0: // nue
      rec_energy = enu - 15.4;
      break;
    case 1:
      rec_energy = enu - 11.4;
      break;
    default:
      std::cerr <<"This sould not be appeared: " << __FILE__ << " L:" << __LINE__ << std::endl;
      rec_energy = 0.0;
      break;
  }

  if(rec_energy > Me){
    totcso = (((crs.at(inifin).at(ienebin)-crs.at(inifin).at(ienebin-1))/(Enu.at(ini).at(ienebin)-Enu.at(ini).at(ienebin-1)))*enu + crs.at(inifin).at(ienebin-1) - ((crs.at(inifin).at(ienebin)-crs.at(inifin).at(ienebin-1))/(Enu.at(ini).at(ienebin)-Enu.at(ini).at(ienebin-1)))*Enu.at(ini).at(ienebin-1)) /* 1.0e-26 */;
  }
  else
    totcso = 0.;
  return std::max(0.0, totcso);
}

std::pair<double,double> SKSNSimXSecNuOxygenSub::GetDiffCrosssection(double e, double r) const 
{
  std::cerr << "Not supported: " << __FILE__ << " GetDiffCrosssection(...) for NuOxygen XSec model (sub)" << std::endl;
  return std::make_pair(0.0,0.0);
}
