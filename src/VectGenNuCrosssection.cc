/**
 * @file VectGenNuCrosssection.cc
 *
 * @date 2017-12-06
 * @author Y.Koshio
2020/05/14: Strumia-Vissani model(2003) was added  BY  M.Harada
 */
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include <math.h>

#include "VectGenSetBin.hh"
#include "VectGenSnConst.hh"
#include "VectGenNuCrosssection.hh"

#include "VectGenUtil.hh"

using namespace std;


VectGenNuCrosssection::VectGenNuCrosssection(){}

static const double costheta_c    = 0.974; // cabibo angle
static const double Delta_R_inner = 0.024; // inner radiative correction
static const double Gf            = 0.04541638E-5; // fermi constant(fm**2), G_F/(hbarc)^3 * (hbarc)^2 in PDG
static const double hbarc         = 197.3269804; // (MeV fm)
static const double f     = 1.;
static const double g     = 1.26; // axial-vector coupling constant


void VectGenNuCrosssection::DcsNuebP_VB(double enu, double costheta, double& Epo, double& dcs)
{
  /*
    get differential cross section of nue_bar P and positron energy
         enu : Neutrino energy in MeV
    costheta : 
         Epo : Positron total energy in MeV
         dcs : Differenctial cross section d sigma / d cos(theta) (cm^2)

    "The angular distribution of the reaction nu_e_bar + p -> e+ + n"
     P.Vogel and J.F.Beacom [arXiv:hep-ph/9903554 1 Apr 1999]
  */

  double Epo0 = enu - DeltaM;

  if(Epo0 < Me) {
    Epo = 0.;
    dcs = 0.;
    return;
  }
  double Pe0 = sqrt(pow(Epo0, 2.) - pow(Me, 2.)); // positron momemtum

  // formula(9)
  double sigma0 = pow(Gf,2.) * pow(costheta_c,2.) / pi * (1. + Delta_R_inner) / pow(hbarc,2.) * 1.0E-26;

  // formula(13)
  double ave_nucleon_mass = (Mp + Mn)/2.; // average nucleon mass
  double v_e_0    = Pe0 / Epo0; // velocity
  double y        = sqrt((pow(DeltaM,2.) - pow(Me,2.))/2.);

  Epo = Epo0 * (1.- enu / ave_nucleon_mass * (1.- v_e_0 *costheta)) - pow(y,2.) / ave_nucleon_mass;

  if(Epo < Me) {
    Epo = 0.;
    dcs = 0.;
    return;
  }

  // formula(15)
  double f2   = 3.706; // nu_p - nu_n (magnetic moment)
  double gamma = 2*(f + f2) * g * ((2*Epo0 + DeltaM) * (1. - v_e_0 * costheta) - pow(Me,2.) / Epo0) + (pow(f,2.) + pow(g,2.)) * (DeltaM * (1. + v_e_0 * costheta) + pow(Me,2.) / Epo0) + (pow(f,2.) + 3 * pow(g,2.)) * ((Epo0 + DeltaM)*(1. - 1./v_e_0 * costheta) - DeltaM) + (pow(f,2.) - pow(g,2.)) * ((Epo0 + DeltaM) * (1. - 1./v_e_0 * costheta) - DeltaM) * v_e_0 * costheta;

  // formula(14)
  double Pe1 = sqrt(pow(Epo, 2.) - pow(Me, 2.)); // positron momemtum
  double v_e_1 = Pe1 / Epo; // velocity

  dcs = sigma0 / 2. * ((pow(f,2.) + 3 * pow(g,2.)) + (pow(f,2.) - pow(g,2.)) * v_e_1 * costheta) * Epo * Pe1 - sigma0 / 2 *(gamma / ave_nucleon_mass) * Epo0 * Pe0;
  //cout << dcs << " " << sigma0 << " " << gamma << " " << v_e_1 << " " << costheta << endl;

  return;
}


void VectGenNuCrosssection::DcsNuebP_SV(double enu, double costheta, double &Epo, double &dcs){ 
  /*********************************************************************/
  /// (purpose)
  ///     Calculate IBD cross section
  /// (input)
  ///     ENU, COSTHETA
  /// (output)
  ///     SIGM, EP
  /// 2019.04 production by H.Ito
  /// 2019.05 inport to slmkb8ibd.F by H.Ito
  /// 2020.05 rewrite to C++ by H.MHarada
  /*********************************************************************/

  double Delta_CM, HBAR_C2, ANO_NU, GF_GeV, GF, G_COUPLING, ALPHA, RHO_NC, totcs, totep;
  //double IBD_THR, ETHRE, IBD_SIG0;

  //     Computations using masses   
  double ave_nucleon_mass = (Mn + Mp)/2.;
  Delta_CM = (Mn*Mn - Mp*Mp - Me*Me)/(2.*Mp);

  //      Fermi coupling constant
  HBAR_C2 = 0.3893793656e-21;// (MeV^2 * cm^2) ... ( plank-constant * photon-velocity )^2
  GF_GeV = 1.1663787e-11;// (MeV^-2)
  GF = GF_GeV*GF_GeV * HBAR_C2;// *1e20*1e20;//(MeV^-2 * cm^2)

  ANO_NU = (2.7928473446-1.)+1.9130427;// sentence btw formula (7) and (8)
  ALPHA = 1. / 137.035999139;

  //MESH_ES = 100;
  //RHO_NC = 1.0126;
  totcs = 0.;
  totep = 0.;


  //IBD_THR  = ((N_MASS+E_MASS)*(N_MASS+E_MASS)- P_MASS*P_MASS)/ (2.*P_MASS);
  //ETHRE = 0.;

  // Positron energy:                                                                                                                                                                    
  double epsilon = enu/Mp;// formula (8)
  double kappa = pow(1.+epsilon,2) - pow(epsilon*costheta,2);// in sentence below formula (21)
  Epo = ((enu-Delta_CM)*(1.+epsilon)+epsilon*costheta*sqrt(pow(enu-Delta_CM,2)-Me*Me*kappa))/kappa;// formula (21)

  // Parameters
  double Ppo = sqrt(Epo*Epo-Me*Me);// formula (21)
  double dE_dCosT = Ppo*epsilon/(1.+epsilon*(1.-costheta*Epo/Ppo));// trans. formula (20)

  // Belows are written in formula (3) - (11)
  double s_minus_u = 2.*Mp*(enu+Epo)-Me*Me;// above of formula (11)
  double t = Mn*Mn-Mp*Mp-2.*Mp*(enu-Epo);// above of formula (11)

  // formula (7)      
  double x  = t /(4.*pow(ave_nucleon_mass,2.)); // part of numerator of 1st formula (7) 
  double y  = 1. - t/710000.;// right side of denominator of 1st formula (7) 
  double z  = 1. - t/1030000.;// denominator of 2nd formula (7)
  //y  = 1. - t/710**2;// right side of denominator of 1st formula (7) 
  //z  = 1. - t/1030**2;// denominator of 2nd formula (7)
  double f1 = (1. - (1.+ANO_NU)*x)/((1.-x)*(y*y));// 1st formula (7)
  double f2 = ANO_NU/((1.-x)*(y*y));// 1st formula (7)
  double g1 = -1.27/z*z;// 2nd formula (7)
  double g2 = 2.*g1*pow(ave_nucleon_mass,2.)/(Mpi*Mpi-t);

  //bolow A,B,C calculation is formula (6)      
  double A = 1./16.*((t-Me*Me)
            *(4.*pow(f1,2.)*(4.*pow(ave_nucleon_mass,2.)+t+Me*Me)
                + 4.*pow(g1,2.)*(-4*pow(ave_nucleon_mass,2.)+t+Me*Me)
                + pow(f2,2.)*(t*t/pow(ave_nucleon_mass,2.)+4.*t+4.*Me*Me)
                + 4*Me*Me*t*pow(g2,2.)/pow(ave_nucleon_mass,2.)
                + 8.*f1*f2*(2*t+Me*Me)
                + 16*Me*Me*g1*g2)
            - pow(DeltaM,2.)
            *((4.*pow(f1,2.)+t*pow(f2,2.)/pow(ave_nucleon_mass,2.))*(4.*pow(ave_nucleon_mass,2.)+t-Me*Me)
                + 4.*pow(g1,2.)*(4.*pow(ave_nucleon_mass,2.)-t+Me*Me)
                + 4.*Me*Me*pow(g2,2.)*(t-Me*Me)/pow(ave_nucleon_mass,2.)
                + 8.*f1*f2*(2.*t-Me*Me)
                + 16.*Me*Me*g1*g2)
            - 32.*Me*Me*ave_nucleon_mass*DeltaM*g1*(f1+f2));

  double B = 1./16.*(16.*t*g1*(f1+f2) + 4.*Me*Me*DeltaM*(pow(f2,2.)+f1*f2+2.*g1*g2)/ave_nucleon_mass);

  double C = 1./16.*(4.*(pow(f1,2)+pow(g1,2)) - t*pow(f2,2.)/pow(ave_nucleon_mass,2.));

  double abs_M_squared = A - B*s_minus_u + C*pow(s_minus_u,2.);// formula (5)
  double rad_correction = ALPHA/pi*(6.00352 + 3./2.*log(Mp/(2.*Epo)) + 1.2*pow(Me/Epo,1.5));// formula (14)        

  double dsigma_dt = GF*pow(costheta_c,2.)/(2*pi*pow(2.*Mp*enu,2.))*abs_M_squared;// formula (3) using s-Mp^2 = 2Mp*Enu(above of formula (11))
  double dsigma_dE = 2.*Mp*dsigma_dt;//formula (11)
  dcs = dE_dCosT*dsigma_dE*(1.+rad_correction);// result 
  //IBD_SIG0 = GF * COSTHETA_C**2 /FUNC_PI * (2. * P_MASS) / (8. * P_MASS**2);
  //SIGM = dE_dCosT * IBD_SIG0 / ENU**2 * abs_M_squared * (1. + rad_correction)  

  if ( dcs < 0 ) dcs = 0.;
  return;
}

double VectGenNuCrosssection::CsNuebP_VB(double enu)
{
  /*
    Total cross section of nu_e_bar + p --> e^+  + n  interaction
    using the precise calculation by Vogel and Beacom
        enu : Neutrino energy in MeV
	(Epo : Positoron total energy in MeV, should be later)
     return : Total cross section (cm^2)
  */

  double totcsnuebp_vb = 0;
  double Ee = enu - DeltaM;
  //if(Ee > Me){
  if(Ee >= 3.){
    for(int j=-100; j<100; j++){
      double costheta = (double(j)+0.5)/100.;
      double dcs0, Epo0;
      VectGenNuCrosssection::DcsNuebP_VB(enu, costheta, Epo0, dcs0);
      totcsnuebp_vb += dcs0;
      //cout << j << " " << Epo0 << " " << dcs0 << " " << totcsnuebp_vb<< endl;
    }
    totcsnuebp_vb /= 100.;
  }

  //cout << totcsnuebp_vb << endl;
  return totcsnuebp_vb;
}

double VectGenNuCrosssection::CsNuebP_SV(double enu)
{
  /*
    Total cross section of nu_e_bar + p --> e^+  + n  interaction
    using the precise calculation by Vogel and Beacom
        enu : Neutrino energy in MeV
	(Epo : Positoron total energy in MeV, should be later)
     return : Total cross section (cm^2)
  */

  double totcsnuebp_SV = 0;
  double Ee = enu - DeltaM;
  if(Ee > Me){
    for(int j=-100; j<100; j++){
      double costheta = (double(j)+0.5)/100.;
      double dcs0, Epo0;
      VectGenNuCrosssection::DcsNuebP_SV(enu, costheta, Epo0, dcs0);
      totcsnuebp_SV += dcs0;
      // cout << j << " " << Epo0 << endl;
    }
    totcsnuebp_SV /= 100.;
  }

  //cout << totcsnuebp_SV << endl;
  return totcsnuebp_SV;
}

void VectGenNuCrosssection::ReadCsNuElastic()
{
  //Open file of total cross-section of elastic scattering
  std::string fnameEla = "/usr/local/sklib_gcc8/skofl-trunk/const/lowe/sn_elastic.root";
  fCsElaFile = new TFile(fnameEla.c_str(), "read");
  fCsElaTree = (TTree*)fCsElaFile->Get("nuela");

  fCsElaTree->SetBranchAddress("nuene", &NuElaEnergy);
  fCsElaTree->SetBranchAddress("cnue", &CsElaNue);
  fCsElaTree->SetBranchAddress("cnux", &CsElaNux);
  fCsElaTree->SetBranchAddress("cneb", &CsElaNeb);
  fCsElaTree->SetBranchAddress("cnxb", &CsElaNxb);
}

double VectGenNuCrosssection::CsNuElastic(int ipart, double enu, int flag)
{
  /*
    Total cross section of nu + e --> nu + e  interaction
    using the precise calculation by Bahcall refer from sollib/sl_elctot_rad.F
        enu : Neutrino energy in MeV
     return : Total cross section (cm^2)
  */

  if((ipart != 12) && (ipart != -12) &&(ipart != 14) && (ipart != -14)){
    std::cerr << "Not support particle code " << ipart << std::endl;
    exit(1);
  }

  double x = 0.;
  if(enu <= nuElaEneMin || enu >= nuElaEneMax) return x;

  if(flag == 0) { // should apply Eth
    double t_min = eEneThr - Me;
    double t_max = 2.*enu*enu/(Me+2.*enu);

    if(t_min > t_max){
      x = 0.;
      return x;
    }

    int istep = 1000;
    double dstep = (t_max - t_min) / double(istep);
    x=0.;
    for (int i = 0; i<istep; i++){
      double E = (t_min + dstep/2.) + dstep * double(i) + Me;
      if(ipart ==  12) x += sl_nue_dif_rad_(&enu, &E) * dstep;
      if(ipart == -12) x += sl_neb_dif_rad_(&enu, &E) * dstep;
      if(ipart ==  14) x += sl_num_dif_rad_(&enu, &E) * dstep;
      if(ipart == -14) x += sl_nmb_dif_rad_(&enu, &E) * dstep;
    }
  }
  else if(flag == 1) { // should not apply Eth
    double e1, e2, c1=0., c2=0.;

    int iene = (int)(enu / nuElaEneBinSize);
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

double VectGenNuCrosssection::calcElectronTotEnergyElastic( const double nuEne, const double cost )
{
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

double VectGenNuCrosssection::calcDeEneDcostElastic( const double nuEne, const double cost )
{
	double alpha = Me / nuEne;

	double p = nuEne * 4. * alpha * cost * SQ( 1 + alpha )
		/ SQ( SQ( 1. + alpha ) - SQ( cost ) );

	return p;
}

double VectGenNuCrosssection::calcCosTthElastic( const double nuE, const double eEth )
{
	double alpha = Me / nuE;
	double y = ( eEth - Me ) / nuE;

	double cosTth = sqrt( y / ( y + 2. * alpha ) ) * ( 1. + alpha );

	return cosTth;
}

