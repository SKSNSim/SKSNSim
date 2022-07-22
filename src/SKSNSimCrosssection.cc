/**************************************
 * File: SKSNSimCrosssection.cc
 * Description:
 * Cross section interface
 *************************************/

#include <cmath>
#include <iostream>
#include "SKSNSimCrosssection.hh"
#include "SKSNSimConstant.hh"

using namespace SKSNSimPhysConst;

std::pair<double,double> SKSNSimXSecIBDSV::GetDiffCrosssection(double enu /* MeV */ , double costheta){
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

  double Delta_CM, G_COUPLING, RHO_NC, totcs, totep;
  //double IBD_THR, ETHRE, IBD_SIG0;

  //     Computations using masses   
  double ave_nucleon_mass = (Mn + Mp)/2.;
  Delta_CM = (Mn*Mn - Mp*Mp - Me*Me)/(2.*Mp);

  //      Fermi coupling constant
  constexpr double HBAR_C2 = 0.3893793656e-21;// (MeV^2 * cm^2) ... ( plank-constant * photon-velocity )^2
  constexpr double GF_GeV = 1.1663787e-11;// (MeV^-2)
  constexpr double GF = GF_GeV*GF_GeV * HBAR_C2;// *1e20*1e20;//(MeV^-2 * cm^2)

  constexpr double ANO_NU = (2.7928473446-1.)+1.9130427;// sentence btw formula (7) and (8)
  constexpr double ALPHA = 1. / 137.035999139;

  //MESH_ES = 100;
  //RHO_NC = 1.0126;
  totcs = 0.;
  totep = 0.;


  //IBD_THR  = ((N_MASS+E_MASS)*(N_MASS+E_MASS)- P_MASS*P_MASS)/ (2.*P_MASS);
  //ETHRE = 0.;

  // Positron energy:                                                                                                                                                                    
  const double epsilon = enu/Mp;// formula (8)
  const double kappa = std::pow(1.+epsilon,2) - std::pow(epsilon*costheta,2);// in sentence below formula (21)
  const double Epo = ((enu-Delta_CM)*(1.+epsilon)+epsilon*costheta*sqrt(pow(enu-Delta_CM,2)-Me*Me*kappa))/kappa;// formula (21)

  // Parameters
  const double Ppo = sqrt(Epo*Epo-Me*Me);// formula (21)
  const double dE_dCosT = Ppo*epsilon/(1.+epsilon*(1.-costheta*Epo/Ppo));// trans. formula (20)

  // Belows are written in formula (3) - (11)
  const double s_minus_u = 2.*Mp*(enu+Epo)-Me*Me;// above of formula (11)
  const double t = Mn*Mn-Mp*Mp-2.*Mp*(enu-Epo);// above of formula (11)

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
  const double A = 1./16.*((t-Me*Me)
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

  const double B = 1./16.*(16.*t*g1*(f1+f2) + 4.*Me*Me*DeltaM*(pow(f2,2.)+f1*f2+2.*g1*g2)/ave_nucleon_mass);

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
