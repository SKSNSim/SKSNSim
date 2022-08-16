/**************************************
 * File: SKSNSimCrosssection.hh
 * Description:
 * Cross section interface
 *************************************/

#ifndef __SKSNSIMCROSSSECTION_H_INCLUDED__
#define __SKSNSIMCROSSSECTION_H_INCLUDED__

#include <utility>
#include <memory>
#include <pdg_codes.h>
#include <TFile.h>
#include <TTree.h>


extern "C" {
	double sl_nue_dif_rad_( double *, double * );
	double sl_neb_dif_rad_( double *, double * );
	double sl_num_dif_rad_( double *, double * );
	double sl_nmb_dif_rad_( double *, double * );
} // TODO to avoid dependency of fortran library

class SKSNSimCrosssectionModel {
  public:
    virtual ~SKSNSimCrosssectionModel() {}
    virtual double /* cm^2 */                          GetCrosssection(double /* MeV */) = 0; // energy -> xsec
    virtual std::pair<double,double> /* <cm^2, MeV> */ GetDiffCrosssection(double /* MeV */, double /* a.u. */) = 0; // energy -> angle -> (xsec, scattered energy)
};

class SKSNSimXSecIBDVB : SKSNSimCrosssectionModel {
  // Cross section model of IBD by Vogel and Beacom
  public:
    SKSNSimXSecIBDVB(){}
    ~SKSNSimXSecIBDVB(){}
    double GetCrosssection(double e);
    std::pair<double,double> GetDiffCrosssection(double e, double r);
};

class SKSNSimXSecIBDSV : SKSNSimCrosssectionModel {
  // Cross section model of IBD by Strumia-Vissani
  public:
    SKSNSimXSecIBDSV(){}
    ~SKSNSimXSecIBDSV(){}
    double GetCrosssection(double);
    std::pair<double,double> GetDiffCrosssection(double, double);
};

class SKSNSimXSecNuElastic : SKSNSimCrosssectionModel {
  // Cross section model of neutrino elastic scattering
  private:
    constexpr static double nuElaEneMin = 0.0; //MeV, minimum neutrino energy 
    constexpr static double nuElaEneMax = 150.0; //MeV, maximum neutrino energy 
    constexpr static int nuElaEneNBins = 15000; // number of bins for nu energy
    constexpr static double nuElaEneBinSize = ( nuElaEneMax - nuElaEneMin ) / ( double )nuElaEneNBins;

    constexpr static double eEneThr = 5.0;// MeV, electron total energy threshold

    const std::string CSELAFILENAME {"/usr/local/sklib_gcc8/skofl-trunk/const/lowe/sn_elastic.root"};
    std::unique_ptr<TFile> fCsElaFile;
    std::unique_ptr<TTree> fCsElaTree;
    double NuElaEnergy, CsElaNue, CsElaNux, CsElaNeb, CsElaNxb;

    void OpenCsElaFile();
    void CloseCsElaFile();
                                          //

  public:
    enum FLAGETHR { ETHRON, ETHROFF };
    SKSNSimXSecNuElastic(){ OpenCsElaFile();}
    ~SKSNSimXSecNuElastic(){CloseCsElaFile();}
    double GetCrosssection(double e, int pid = -PDG_ELECTRON_NEUTRINO, FLAGETHR flag = ETHRON);
    std::pair<double,double> GetDiffCrosssection(double e, double r);
    static double GetNuEneMin() {return nuElaEneMin;}
    static double GetNuEneMax() {return nuElaEneMax;}
};

class SKSNSimXSecNuOxygen : SKSNSimCrosssectionModel {
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
    std::map<INISTATE,int> exNum; // = {{3,15,8,1,6}};
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

    void LoadFile(INISTATE);
    void InitializeTable();

  public:
    SKSNSimXSecNuOxygen(){
      InitializeTable();
    }
    ~SKSNSimXSecNuOxygen(){}
    double GetCrosssection(double e, INIFINSTATE inifin = {0,0,0,0});
    std::pair<double,double> GetDiffCrosssection(double, double);
};

class SKSNSimXSecNuOxygenSub : SKSNSimCrosssectionModel {
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

    void LoadFile(INISTATE);
    void InitializeTable();

  public:
    SKSNSimXSecNuOxygenSub(){
      InitializeTable();
    }
    ~SKSNSimXSecNuOxygenSub(){}
    double GetCrosssection(double e, INIFINSTATE inifin = {0,0,0});
    std::pair<double,double> GetDiffCrosssection(double, double);
};

#endif
