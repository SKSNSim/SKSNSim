/*********************************
 * SKSNSimConstant.hh
 ********************************/
#include <cmath>
#include <map>
#include "SKSNSimEnum.hh"

namespace SKSNSimPhysConst {
  constexpr double Me = 0.510998950e0 /* MeV */;// electron mass(from PDG 2020) // TODO this is different against with Me in sl_nue_dif_rad.F in sollib
  constexpr double Mp = 938.272088e0 /* MeV */; // proton mass
  constexpr double Mn = 939.565421e0 /* MeV */; // neutron mass
  constexpr double Mpi = 139.57061e0 /* MeV */; // pion mass
  constexpr double DeltaM = Mn - Mp; //mass difference proton and neutron
  constexpr double costheta_cabibo = 0.974; // cabibo angle
  constexpr double Gf = 0.04541638e-5; // fermi constant(fm**2), G_F/(hbarc)^3 * (hbarc)^2 in PDG
  constexpr double HBARC = 197.3269804; // (MeV fm)

  constexpr double PI = M_PI /*3.14159265358979833*/;
  constexpr double ERG2MEV = 1.0e0/1.60217733e-6;

  constexpr double ZERO_PRECISION = 1e-9;

  // Number of target in 32.48 kton of SK total Volume
  constexpr double Ntarget_p = 2.173e33;
  constexpr double Ntarget_e = 1.086e34;
  constexpr double Ntarget_o = 1.086e33;

  // Unit of the distance to the Supenova 10kpc = 3.086 x 10^22 cm
  constexpr double Distance = 3.086e22; // cm for 10kpc

  constexpr double Const_p = 1.0 / (4.0*M_PI*Distance*Distance) * Ntarget_p;
  constexpr double Const_e = 1.0 / (4.0*M_PI*Distance*Distance) * Ntarget_e;
  constexpr double Const_o = 1.0 / (4.0*M_PI*Distance*Distance) * Ntarget_o;


  //========================
  // Neutrino Oscillation Parameter
  constexpr double PMNSSinSqTheta12 = 0.28; // PMNS matrix
  constexpr double PMNSCosSqTheta12 = 1. - PMNSSinSqTheta12;
  const std::map<SKSNSIMENUM::NEUTRINOOSCILLATION, std::tuple<double,double,double,double,double,double,double,double>> NuOscProbCollection = {
    /* Oscillation type -> { nue1, nue2, neb1, neb2, nux1, nux2, nxb1, nxb2} */
    {SKSNSIMENUM::NEUTRINOOSCILLATION::kNONE,     { 1., 0., 1., 0., 2., 0., 2., 0.}},
    {SKSNSIMENUM::NEUTRINOOSCILLATION::kNORMAL,   { 0., 1., PMNSCosSqTheta12, PMNSSinSqTheta12, 1., 1., 1.+PMNSCosSqTheta12, PMNSSinSqTheta12}},
    {SKSNSIMENUM::NEUTRINOOSCILLATION::kINVERTED, { PMNSSinSqTheta12, PMNSCosSqTheta12, 0., 1., 1.+PMNSSinSqTheta12, PMNSCosSqTheta12, 1., 1.}}
  };
  // Acsessing Functions
  double GetNuOscNue1(SKSNSIMENUM::NEUTRINOOSCILLATION t){ return std::get<0>(NuOscProbCollection.at(t));}
  double GetNuOscNue2(SKSNSIMENUM::NEUTRINOOSCILLATION t){ return std::get<1>(NuOscProbCollection.at(t));}
  double GetNuOscNueb1(SKSNSIMENUM::NEUTRINOOSCILLATION t){return std::get<2>(NuOscProbCollection.at(t));}
  double GetNuOscNueb2(SKSNSIMENUM::NEUTRINOOSCILLATION t){return std::get<3>(NuOscProbCollection.at(t));}
  double GetNuOscNux1(SKSNSIMENUM::NEUTRINOOSCILLATION t){ return std::get<4>(NuOscProbCollection.at(t));}
  double GetNuOscNux2(SKSNSIMENUM::NEUTRINOOSCILLATION t){ return std::get<5>(NuOscProbCollection.at(t));}
  double GetNuOscNuxb1(SKSNSIMENUM::NEUTRINOOSCILLATION t){return std::get<6>(NuOscProbCollection.at(t));}
  double GetNuOscNuxb2(SKSNSIMENUM::NEUTRINOOSCILLATION t){return std::get<7>(NuOscProbCollection.at(t));}

}
