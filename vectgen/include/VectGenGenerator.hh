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

#include "VectGenSnRoot.hh"

/**
 * @class Generator
 */
class VectGenGenerator
{
public:
  VectGenGenerator(){}
  ~VectGenGenerator(){}

  void convDirection(const double, const double, double*);
  void determineAngleNuebarP(const double, double&, double&, double&);
  void determineAngleElastic(const int, const double, double&, double&, double&);
  void determineKinematics(const int, const double, double*, MCInfo*);
  void determinePosition(double&, double&, double&);
  void FillEvent();
  void MakeEvent(double, double, int, int, double);
  void Process();
  void Process(int);

protected:

  TRandom3 *generator;
  int flag_event;
  double RatioTo10kpc;
  std::string OutDir;
  double snDir[3];
  double Rmat[3][3];

  VectGenSnFlux* nuflux;
  VectGenNuCrosssection* nucrs;

  //Neutrino oscillation
  double oscnue1, oscnue2, oscneb1, oscneb2, oscnux1, oscnux2, oscnxb1, oscnxb2;

  //expected total number of events
  double totNuebarp;
  double totNueElastic, totNuebarElastic, totNuxElastic, totNuxbarElastic;
  double totNueO, totNuebarO;
  double totNcNup, totNcNun, totNcNubarp, totNcNubarn;

private:
  //store neutrino kinematics into vector for time sorting
  int nReact, nuType;
  std::vector<SNEvtInfo> vEvtInfo;

  // number of events where file switches
  const int NeventFile = 50000;

  int iSkip;
};

#endif
