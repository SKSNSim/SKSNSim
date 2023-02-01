#ifndef VECTGENSNCONST_H_INCLUDED
#define VECTGENSNCONST_H_INCLUDED

// basic parameters
const double Me = 0.510998950e0;// electron mass(from PDG 2020)
const double Mp = 938.272088e0; // proton mass
const double Mn = 939.565421e0; // neutron mass
const double Mpi = 139.57039; // pion mass
 
//const double DeltaM = 1.2933317e0; //mass difference proton and neutron
double DeltaM = Mn - Mp; //mass difference proton and neutron
const double pi = 3.14159265358979833;

// neutrino oscillation
#ifdef ORIGINAL_NUOSCPARAMETER
const double sin2th12  = 0.28;
#else
const double sin2th12  = 0.307; // From PDG 2022 (R.L. Workmanet al.(Particle Data Group), Prog.Theor.Exp.Phys.2022, 083C01 (2022))
#endif
double cos2th12 = 1. - sin2th12;

//
// Unit is 32.48 kton in 10kpc for making time-spectrum table
//

// Number of target in 32.48 kton of SK total Volume
const double Ntarget_p = 2.173e33;
const double Ntarget_e = 1.086e34;
const double Ntarget_o = 1.086e33;

// Unit of the distance to the Supenova 10kpc = 3.086 x 10^22 cm
const double Distance = 3.086e22; // cm for 10kpc
//const double DistanceUnit = 3.086e18; // cm for 1pc

double Const_p = 1.0 / (4.0*pi*pow(Distance,2.)) * Ntarget_p;
double Const_e = 1.0 / (4.0*pi*pow(Distance,2.)) * Ntarget_e;
double Const_o = 1.0 / (4.0*pi*pow(Distance,2.)) * Ntarget_o;

#endif // SN_CONST_H_INCLUDED
