/**************************************
 * File: SKSNSimCrosssection.hh
 * Description:
 * Cross section interface
 *************************************/

#ifndef SKSNSIMCROSSSECTION_H_INCLUDED
#define SKSNSIMCROSSSECTION_H_INCLUDED

#include <utility>
#include <memory>
#include <set>
#include <pdg_codes.h>
#include <TFile.h>
#include <TTree.h>
#include "SKSNSimConstant.hh"


extern "C" {
	double sl_nue_dif_rad_( double *, double * );
	double sl_neb_dif_rad_( double *, double * );
	double sl_num_dif_rad_( double *, double * );
	double sl_nmb_dif_rad_( double *, double * );
} // TODO to avoid dependency of fortran library

enum struct XSECTYPE { mXSECIBD = 0, mXSECELASTIC, mXSECOXYGEN, mXSECOXYGENSUB, mXSECOXYGENNC, mNXSECTYPE};
enum struct NUREACTTYPE { kNUEBARP = 0,
  kNUEELASTIC, kNUEBARELASTIC, kNUXELASTIC, kNUXBARELASTIC,
  kNUEO, kNUEBARO,
  kNUEPNC, kNUENNC, kNUEBARPNC, kNUEBARNNC, kNUXPNC, kNUXNNC, kNUXBARPNC, kNUXBARNNC,
  kNNUREACTTYPE};

namespace SKSNSimCrosssection {
  double CalcIBDEnuFromEpos ( const double /* positron energy MeV */, const double /* cosTheta between positron and neutrino */);
}

class SKSNSimCrosssectionModel {
  public:
    enum XSECNUTYPE { XSECNUE = 0, XSECNUEB, XSECNUX, NXSECNUTYPE};
    virtual ~SKSNSimCrosssectionModel() {}
    virtual double /* cm^2 */                          GetCrosssection(double /* MeV */) const = 0; // energy -> xsec
    virtual std::pair<double,double> /* <cm^2, MeV> */ GetDiffCrosssection(double /* MeV */, double /* a.u. */) const = 0; // energy -> angle -> (xsec, scattered energy)
    //virtual const std::set<XSECNUTYPE> &GetSupportedNuType () const = 0;

};

class SKSNSimXSecFlat : public SKSNSimCrosssectionModel {
  // Flat cross section: always return 1.0
  private:
  public:
    SKSNSimXSecFlat(){}
    ~SKSNSimXSecFlat(){}
    double GetCrosssection(double e){ return 1.0;}
    std::pair<double,double> GetDiffCrosssection(double e, double r) { return std::make_pair(1.0, e);};
};


class SKSNSimXSecIBDVB : public SKSNSimCrosssectionModel {
  // Cross section model of IBD by Vogel and Beacom
  private:
  public:
    SKSNSimXSecIBDVB(){}
    ~SKSNSimXSecIBDVB(){}
    double GetCrosssection(double e) const;
    std::pair<double,double> GetDiffCrosssection(double e, double r) const;
};

class SKSNSimXSecIBDSV : public SKSNSimCrosssectionModel {
  // Cross section model of IBD by Strumia-Vissani
  private:
  public:
    SKSNSimXSecIBDSV(){}
    ~SKSNSimXSecIBDSV(){}
    double GetCrosssection(double) const;
    std::pair<double,double> GetDiffCrosssection(double, double) const;
};

class SKSNSimXSecIBDRVV : public SKSNSimCrosssectionModel {
  // Cross section model of IBD by Ricciardi-Vignaroli-Vissani (DOI: https://doi.org/10.1007/JHEP08(2022)212)
  // Reference [2] IBD calculatoin of Strumia-Vissani (Phys.Lett.B564:42-54,2003, DOI: https://doi.org/10.1016/S0370-2693(03)00616-6)
  // Error is systematic uncertainty cased by (Vud, axial coupling) and axial_radii. This is implemented as simple approximation with constant and power-law, respectively.
  private:
  public:
    SKSNSimXSecIBDRVV(){}
    virtual ~SKSNSimXSecIBDRVV(){}
    double GetCrosssection(double) const;
    double GetCrosssectionError(double) const; /* arbitary unit. If we want to convert to unit of %, multiply 100.0. */
    std::pair<double,double> GetDiffCrosssection(double, double) const;
};

class SKSNSimXSecNuElastic : public SKSNSimCrosssectionModel {
  // Cross section model of neutrino elastic scattering
  private:
    constexpr static double nuElaEneMin = 0.0; //MeV, minimum neutrino energy 
    constexpr static double nuElaEneMax = 150.0; //MeV, maximum neutrino energy 
    constexpr static int nuElaEneNBins = 15000; // number of bins for nu energy
    constexpr static double nuElaEneBinSize = ( nuElaEneMax - nuElaEneMin ) / ( double )nuElaEneNBins;

    constexpr static double eEneThr = 5.0;// MeV, electron total energy threshold

    std::unique_ptr<TFile> fCsElaFile;
    TTree *fCsElaTree;
    double NuElaEnergy, CsElaNue, CsElaNux, CsElaNeb, CsElaNxb;

    void OpenCsElaFile();
    void CloseCsElaFile();

  public:
    enum FLAGETHR { ETHRON, ETHROFF };
    SKSNSimXSecNuElastic(){ OpenCsElaFile();}
    ~SKSNSimXSecNuElastic(){CloseCsElaFile();}
    double GetCrosssection(double e, int pid = -PDG_ELECTRON_NEUTRINO, FLAGETHR flag = ETHRON) const;
    double GetCrosssection(double e) const { return GetCrosssection(e, -PDG_ELECTRON_NEUTRINO, ETHRON);} ;
    std::pair<double,double> GetDiffCrosssection(double e, double r) const;
    static double GetNuEneMin() {return nuElaEneMin;}
    static double GetNuEneMax() {return nuElaEneMax;}
    static double CalcElectronTotEnergy( const double , const double );
    static double CalcDeEneDCost( const double , const double );
    static double CalcCosThr( const double , const double );
};

class SKSNSimXSecNuOxygen : public SKSNSimCrosssectionModel {
  // Cross section model of neutrino-oxygen -> single lepton
  private:
    typedef std::tuple<int,int> INISTATE; // <nue-or-nuebar,ix>
    typedef std::tuple<int,int,int,int> INIFINSTATE; // <nue-or-nuebar,ix,ex,ch>
    constexpr static int NTYPE = 2;
    constexpr static int NIXSTATE = 5;
    constexpr static int NEXSTATE = 16; // maximum, depending to IXSTATE
    constexpr static int NCHANNEL = 7;
    constexpr static double eEneThr = 5.0;// MeV, electron total energy threshold
    std::map<INISTATE,bool> isOpened;
    std::map<INISTATE, std::vector<double>> rec;
    std::map<INISTATE, std::vector<double>> pro_energy;
    std::map<INISTATE,int> fileSize;
    std::map<INISTATE,int> exNum;
    std::map<INISTATE,std::vector<double>> exEne;
    std::map<INISTATE,std::vector<double>> nuene;
    std::map<INIFINSTATE, std::vector<double>> crs;


    static int GetNumEx(int ix){
      const static int num_ex[NIXSTATE] = {3, 15, 8, 1, 16};
      return num_ex[ix];
    }
    static INISTATE    convToINISTATE(int type, int ix) {return std::make_tuple(type,ix);}
    static INIFINSTATE convToINIFINSTATE(int type, int ix, int ex, int ch) {return std::make_tuple(type,ix,ex,ch);}
    static INISTATE GetINISTATE(INIFINSTATE f){ return convToINISTATE(std::get<0>(f), std::get<1>(f));}
    static INIFINSTATE convToINIFINSTATE(INISTATE ini, int ex, int ch) {return std::make_tuple(std::get<0>(ini),std::get<1>(ini),ex,ch);}

    void LoadFile();
    void LoadFile(INISTATE ini);
    void InitializeTable();
  public:
    SKSNSimXSecNuOxygen(){
      InitializeTable();
      LoadFile();
    }
    ~SKSNSimXSecNuOxygen(){}
    double GetCrosssection(double e, INIFINSTATE inifin = {0,0,0,0}) const;
    double GetCrosssection(double e) const { return GetCrosssection(e, {0,0,0,0});};
    std::pair<double,double> GetDiffCrosssection(double, double) const;

    double OxigFuncAngleRecCC(int num, int ix, int ex, int ch, double enu, double cos);
    double OxigFuncRecEneCC(int num, int ix, int ex, int ch, double enu);
};

class SKSNSimXSecNuOxygenNC : public SKSNSimCrosssectionModel {
  // Cross section model of neutrino-oxygen -> single lepton via NC
  private:
    typedef std::tuple<int> INISTATE; // <nue-or-nuebar,ix>
    typedef std::tuple<int,int> INIFINSTATE; // <nue-or-nuebar,ix,ex,ch>
    constexpr static int NTYPE = 2;
    constexpr static int NEXSTATE = 8; // maximum, depending to IXSTATE, TYPE
    // constexpr static double eEneThr = 5.0;// MeV, electron total energy threshold
    std::map<INISTATE,bool> isOpened;
    std::map<INISTATE,int> fileSize;
    std::map<INISTATE,std::vector<double>> exEne;
    std::map<INISTATE,std::vector<double>> nuene;
    std::map<INIFINSTATE, std::vector<double>> crs;


    static INISTATE    convToINISTATE(int type) {return std::make_tuple(type);}
    static INIFINSTATE convToINIFINSTATE(int type, int ex) {return std::make_tuple(type,ex);}
    static INISTATE GetINISTATE(INIFINSTATE f){ return convToINISTATE(std::get<0>(f));}
    static INIFINSTATE convToINIFINSTATE(INISTATE ini, int ex) {return std::make_tuple(std::get<0>(ini),ex);}

    void LoadFile();
    void LoadFile(INISTATE ini);
    void InitializeTable();
  public:
    SKSNSimXSecNuOxygenNC(){
      InitializeTable();
      LoadFile();
    }
    ~SKSNSimXSecNuOxygenNC(){}
    double GetCrosssection(double e, INIFINSTATE inifin = {0,0}) const;
    double GetCrosssection(double e) const { return GetCrosssection(e, {0,0});};
    std::pair<double,double> GetDiffCrosssection(double, double) const;
    static int GetNumEx(int type){
      const static int num_ex[NTYPE] = {8, 4};
      return num_ex[type];
    }

};

class SKSNSimXSecNuOxygenSub : public SKSNSimCrosssectionModel {
  // Cross section model of neutrino-oxygen -> single lepton
  private:
    typedef std::tuple<int,int> INISTATE; // <nue-or-nuebar,ix>
    typedef std::tuple<int,int,int> INIFINSTATE; // <nue-or-nuebar,ix,ch>
    constexpr static int NTYPE = 2;
    constexpr static int NIXSTATE = 5;
    constexpr static int NCHANNEL = 32;
    constexpr static double eEneThr = 5.0;// MeV, electron total energy threshold
    std::map<INISTATE,bool> isOpened;
    std::map<INISTATE,int> fileSize;
    std::map<INISTATE,std::vector<double>> Enu;
    std::map<INIFINSTATE, std::vector<double>> crs;

    static INISTATE    convToINISTATE(int type, int ix) {return std::make_tuple(type,ix);}
    static INIFINSTATE convToINIFINSTATE(int type, int ix, int ch) {return std::make_tuple(type,ix,ch);}
    static INISTATE GetINISTATE(INIFINSTATE f){ return convToINISTATE(std::get<0>(f), std::get<1>(f));}
    static INIFINSTATE convToINIFINSTATE(INISTATE ini, int ch) {return std::make_tuple(std::get<0>(ini),std::get<1>(ini),ch);}

    void LoadFile();
    void LoadFile(INISTATE);
    void InitializeTable();


  public:
    SKSNSimXSecNuOxygenSub(){
      InitializeTable();
      LoadFile();
    }
    ~SKSNSimXSecNuOxygenSub(){}
    double GetCrosssection(double e, INIFINSTATE inifin = {0,0,0}) const;
    double GetCrosssection(double e) const { return GetCrosssection(e, {0,0,0});};
    std::pair<double,double> GetDiffCrosssection(double, double) const;

};

#endif
