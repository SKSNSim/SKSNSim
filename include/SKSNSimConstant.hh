/*********************************
 * SKSNSimConstant.hh
 ********************************/

namespace SKSNSimPhysConst {
constexpr double Me = 0.510998950e0 /* MeV */;// electron mass(from PDG 2020)
constexpr double Mp = 938.272088e0 /* MeV */; // proton mass
constexpr double Mn = 939.565421e0 /* MeV */; // neutron mass
constexpr double Mpi = 139.57061e0 /* MeV */; // pion mass
constexpr double DeltaM = Mn - Mp; //mass difference proton and neutron
constexpr double costheta_cabibo = 0.974; // cabibo angle
constexpr double Gf = 0.04541638E-5; // fermi constant(fm**2), G_F/(hbarc)^3 * (hbarc)^2 in PDG
constexpr double HBARC = 197.3269804; // (MeV fm)

constexpr double PI = 3.14159265358979833;
 
}
