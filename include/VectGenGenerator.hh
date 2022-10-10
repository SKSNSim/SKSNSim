/**
 * @file VectGenGenerator.hh
 *
 * @date 2021-09-22
 * @author Y.Koshio
 */

#ifndef VECTGENGENERATOR_HH
#define VECTGENGENERATOR_HH

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <algorithm>

#include "snevtinfo.h"
#include "mcinfo.h"
#include "geotnkC.h"

#include "VectGenSetBin.hh"
#include "VectGenSnConst.hh"
#include "VectGenSnFlux.hh"
#include "VectGenNuCrosssection.hh"
#include "VectGenOxyCrosssection.hh"
#include "VectGenOxyCrosssection_sub.hh"
#include "VectGenOxigFunc.hh"

#include "VectGenSnRoot.hh"

#include "FluxCalculation.hh"

/**
 * @class Generator
 */
extern "C" {
  void read_timevent_( int *, int *, int*);
}

extern "C" {
	void dir_oxigfunc_( double *, double *, double *, double * );
}
enum PositionType{
  mInnerFV = 0,  // 0
  mInnerID,  // 1
  mEntireTank  // 2
} ;

enum SKRUN
{
    SK_IV_BEGIN = 60000,
    SK_IV_END = 79999,
    SK_V_BEGIN  = 80000,
    SK_V_END  = 84999,
    SK_VI_BEGIN  = 85000,
};

class VectGenGenerator
{
public:
  VectGenGenerator();
  ~VectGenGenerator(){}

  void convDirection(const double, const double, double*);
  void determineAngleNuebarP(const double, double&, double&, double&);
  void determineAngleElastic(const int, const double, double&, double&, double&);
  void determineAngleNueO(const int, const int, const int, const int, const double, double&, double&, double&);
  void determineKinematics(const int, const double, double*, MCInfo*);
  void determinePosition(int, double&, double&, double&);
  void determineNmomentum(double&, double&, double&);
  void FillEvent();
  void MakeEvent(double, double, int, int, double);
  void Process();    // For SN generator
  void Process(int); // For DSBN vector generator
  void ReadTimeEventFile(int* numevent, int subrun[]);

protected:

  TRandom3 *generator;
  int flag_event;
  double RatioTo10kpc;
  std::string OutDir;
  std::string ModelName;
  double snDir[3];
  double Rmat[3][3];

  VectGenSnFlux* nuflux;
  VectGenNuCrosssection* nucrs;
  std::unique_ptr<FluxCalculation> nuflux_dsnb;
  VectGenOxyCrosssection* ocrs;
  VectGenOxyCrosssectionSub* osub;
  VectGenOxyCrosssection* ocrs_nc;
  VectGenOxigFunc* reco;
  VectGenOxigFunc* rece;

  //Neutrino oscillation
  double oscnue1, oscnue2, oscneb1, oscneb2, oscnux1, oscnux2, oscnxb1, oscnxb2;

  //expected total number of events
  double totNuebarp;
  double totNueElastic, totNuebarElastic, totNuxElastic, totNuxbarElastic;
  double totNueO, totNuebarO;
  double totNueOsub, totNuebarOsub;
  double totNcNuep, totNcNuebarp, totNcNuxp, totNcNuxbarp, totNcNuen, totNcNuebarn, totNcNuxn, totNcNuxbarn;
  double totNcNuep0, totNcNuebarp0, totNcNuxp0, totNcNuxbarp0;
  double totNcNuep1, totNcNuebarp1, totNcNuxp1, totNcNuxbarp1;
  double totNcNuep2, totNcNuebarp2, totNcNuxp2, totNcNuxbarp2;
  double totNcNuep3, totNcNuebarp3, totNcNuxp3, totNcNuxbarp3;
  double totNcNuep4, totNcNuebarp4, totNcNuxp4, totNcNuxbarp4;
  double totNcNuep5, totNcNuebarp5, totNcNuxp5, totNcNuxbarp5;
  double totNcNuep6, totNcNuebarp6, totNcNuxp6, totNcNuxbarp6;
  double totNcNuep7, totNcNuebarp7, totNcNuxp7, totNcNuxbarp7;
  double totNcNuen0, totNcNuebarn0, totNcNuxn0, totNcNuxbarn0;
  double totNcNuen1, totNcNuebarn1, totNcNuxn1, totNcNuxbarn1;
  double totNcNuen2, totNcNuebarn2, totNcNuxn2, totNcNuxbarn2;
  double totNcNuen3, totNcNuebarn3, totNcNuxn3, totNcNuxbarn3;

  int totGenNuebarp=0;
  int totGenNueElastic=0, totGenNuebarElastic=0, totGenNuxElastic=0, totGenNuxbarElastic=0;
  int totGenNueO=0, totGenNuebarO=0;
  int totGenNueOsub=0, totGenNuebarOsub=0;
  //int totGenNcNup=0, totGenNcNun=0, totGenNcNubarp=0, totGenNcNubarn=0;
  int totGenNcNuep=0, totGenNcNuebarp=0, totGenNcNuxp=0, totGenNcNuxbarp=0, totGenNcNuen=0, totGenNcNuebarn=0, totGenNcNuxn=0, totGenNcNuxbarn=0;
  int totGenNcNuep0=0, totGenNcNuebarp0=0, totGenNcNuxp0=0, totGenNcNuxbarp0=0;
  int totGenNcNuep1=0, totGenNcNuebarp1=0, totGenNcNuxp1=0, totGenNcNuxbarp1=0;
  int totGenNcNuep2=0, totGenNcNuebarp2=0, totGenNcNuxp2=0, totGenNcNuxbarp2=0;
  int totGenNcNuep3=0, totGenNcNuebarp3=0, totGenNcNuxp3=0, totGenNcNuxbarp3=0;
  int totGenNcNuep4=0, totGenNcNuebarp4=0, totGenNcNuxp4=0, totGenNcNuxbarp4=0;
  int totGenNcNuep5=0, totGenNcNuebarp5=0, totGenNcNuxp5=0, totGenNcNuxbarp5=0;
  int totGenNcNuep6=0, totGenNcNuebarp6=0, totGenNcNuxp6=0, totGenNcNuxbarp6=0;
  int totGenNcNuep7=0, totGenNcNuebarp7=0, totGenNcNuxp7=0, totGenNcNuxbarp7=0;
  int totGenNcNuen0=0, totGenNcNuebarn0=0, totGenNcNuxn0=0, totGenNcNuxbarn0=0;
  int totGenNcNuen1=0, totGenNcNuebarn1=0, totGenNcNuxn1=0, totGenNcNuxbarn1=0;
  int totGenNcNuen2=0, totGenNcNuebarn2=0, totGenNcNuxn2=0, totGenNcNuxbarn2=0;
  int totGenNcNuen3=0, totGenNcNuebarn3=0, totGenNcNuxn3=0, totGenNcNuxbarn3=0;
 

  // parameters
  double nuEne_min, nuEne_max, maxProb;

  //number of particle emitted on deexcitation with CC reaction
  int numNtNueO[7] = {0, 1, 0, 2, 0, 0, 0};
  int numPtNueO[7] = {1, 1, 2, 1, 0, 0, 1};
  int numNtNuebarO[7] = {0, 1, 0, 2, 1, 0, 0};
  int numPtNuebarO[7] = {0, 0, 1, 0, 1, 1, 2};
  int numGmNuebarO[7] = {1, 0, 0, 0, 0, 0, 0};

  //energy of gamma on deexcitation with NC reaction
  double eneGamN[8] = {5.27, 6.33, 7.16, 7.56, 8.32, 8.57, 9.05, 9.76};
  double eneGamO[4] = {5.18, 6.18, 6.69, 7.28};
  //flavor of neutrinos
  int neutrinoType[4] = {12,-12,14,-14};

  std::vector<double> Ocrse0[16][7];
  std::vector<double> Ocrse1[16][7];
  std::vector<double> Ocrse2[16][7];
  std::vector<double> Ocrse3[16][7];
  std::vector<double> Ocrse4[16][7];
  std::vector<double> Ocrsp0[16][7];
  std::vector<double> Ocrsp1[16][7];
  std::vector<double> Ocrsp2[16][7];
  std::vector<double> Ocrsp3[16][7];
  std::vector<double> Ocrsp4[16][7];
  std::vector<double> OcrseSub[5][32];
  std::vector<double> OcrspSub[5][32]; 
  std::vector<double> OcrsNC[2][14];
  //std::vector<double> OcrsNC[2];

  // MC data
  TFile* fOutFile;
  TTree* theOTree;

  int fRefRunNum;
  bool bIsUseTimeEvent;
  bool bUseFlatFlux;

  MCInfo* fMC = 0;

private:
  //store neutrino kinematics into vector for time sorting
  int nReact, nuType;
  std::string sReact;
  std::vector<SNEvtInfo> vEvtInfo;

  // number of events where file switches
  const int NeventFile = 50000;

  int iSkip;
};

#endif
