/**************************************
 * File: SKSNSimCrosssection.hh
 * Description:
 * Cross section interface
 *************************************/

#ifndef __SKSNSIMCROSSSECTION_H_INCLUDED__
#define __SKSNSIMCROSSSECTION_H_INCLUDED__

#include <utility>

class SKSNSimCrosssectionModel {
  public:
    virtual ~SKSNSimCrosssectionModel() {}
    virtual double /* cm^2 */                          GetCrosssection(double /* MeV */) = 0; // energy -> xsec
    virtual std::pair<double,double> /* <cm^2, MeV> */ GetDiffCrosssection(double /* MeV */, double /* a.u. */) = 0; // energy -> angle -> (xsec, scattered energy)
};

class SKSNSimXSecIBDVB : SKSNSimCrosssectionModel {
  // Cross section model of IBD by VB
  public:
    SKSNSimXSecIBDVB(){}
    ~SKSNSimXSecIBDVB(){}
    double GetCrosssection(double e) { return 0.0;}
    std::pair<double,double> GetDiffCrosssection(double e, double r) { return std::make_pair(0.0,0.0);}
};

class SKSNSimXSecIBDSV : SKSNSimCrosssectionModel {
  // Cross section model of IBD by Strumia-Vissani
  public:
    SKSNSimXSecIBDSV(){}
    ~SKSNSimXSecIBDSV(){}
    double GetCrosssection(double e) { return 0.0;}
    std::pair<double,double> GetDiffCrosssection(double, double);
};

class SKSNSimXSecNuElastic : SKSNSimCrosssectionModel {
  // Cross section model of neutrino elastic scattering
  public:
    SKSNSimXSecNuElastic(){}
    ~SKSNSimXSecNuElastic(){}
    double GetCrosssection(double e) { return 0.0;}
    std::pair<double,double> GetDiffCrosssection(double e, double r) { return std::make_pair(0.0,0.0);}
};

class SKSNSimXSecNuOxygen : SKSNSimCrosssectionModel {
  // Cross section model of neutrino-oxygen -> single lepton
  public:
    SKSNSimXSecNuOxygen(){}
    ~SKSNSimXSecNuOxygen(){}
    double GetCrosssection(double e) { return 0.0;}
    std::pair<double,double> GetDiffCrosssection(double e, double r) { return std::make_pair(0.0,0.0);}
};

#endif
