/**************************************
 * File: SKSNSimCrosssection.hh
 * Description:
 * Cross section interface
 *************************************/

#ifndef __SKSNSIMCROSSSECTION_H_INCLUDED__
#define __SKSNSIMCROSSSECTION_H_INCLUDED__

#include <utility>
#include <memory>
#include <set>
#include <pdg_codes.h>
#include <TFile.h>
#include <TTree.h>


extern "C" {
	double sl_nue_dif_rad_( double *, double * );
	double sl_neb_dif_rad_( double *, double * );
	double sl_num_dif_rad_( double *, double * );
	double sl_nmb_dif_rad_( double *, double * );
} // TODO to avoid dependency of fortran library

enum XSECTYPE { mXSECIBD = 0, mXSECELASTIC, mXSECOXYGEN, mXSECOXYGENSUB, mNXSECTYPE};

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

class SKSNSimXSecNuElastic : public SKSNSimCrosssectionModel {
  // Cross section model of neutrino elastic scattering
  private:
    constexpr static double nuElaEneMin = 0.0; //MeV, minimum neutrino energy 
    constexpr static double nuElaEneMax = 150.0; //MeV, maximum neutrino energy 
    constexpr static int nuElaEneNBins = 15000; // number of bins for nu energy
    constexpr static double nuElaEneBinSize = ( nuElaEneMax - nuElaEneMin ) / ( double )nuElaEneNBins;

    constexpr static double eEneThr = 5.0;// MeV, electron total energy threshold

    const std::string CSELAFILENAME {"/usr/local/sklib_gcc8/skofl-trunk/const/lowe/sn_elastic.root"};
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
    typedef std::tuple<int,int> INISTATE; // <nue-or-nuebar,ix>
    typedef std::tuple<int,int,int,int> INIFINSTATE; // <nue-or-nuebar,ix,ex,ch>
    constexpr static int NTYPE = 2;
    constexpr static int NIXSTATE = 2;
    constexpr static int NEXSTATE = 8; // maximum, depending to IXSTATE
    constexpr static int NCHANNEL = 1;
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
      const static int num_ex[NIXSTATE] = {8, 4};
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
    SKSNSimXSecNuOxygenNC(){
      InitializeTable();
      LoadFile();
    }
    ~SKSNSimXSecNuOxygenNC(){}
    double GetCrosssection(double e, INIFINSTATE inifin = {0,0,0,0}) const;
    double GetCrosssection(double e) const { return GetCrosssection(e, {0,0,0,0});};
    std::pair<double,double> GetDiffCrosssection(double, double) const;

    double OxigFuncAngleRecCC(int num, int ix, int ex, int ch, double enu, double cos);
    double OxigFuncRecEneCC(int num, int ix, int ex, int ch, double enu);
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
