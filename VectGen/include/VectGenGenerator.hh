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

  void Process();

protected:

  VectGenSnFlux* nuflux;
  VectGenNuCrosssection* nucrs;

  TH1D* evrate;

  //Neutrino oscillation
  double oscnue1, oscnue2, oscneb1, oscneb2, oscnux1, oscnux2, oscnxb1, oscnxb2;

  //expected total number of events
  double totNuebarp;
  double totNueElastic, totNuebarElastic, totNuxElastic, totNuxbarElastic;
  double totNueO, totNuebarO;
  double totNcNup, totNcNun, totNcNubarp, totNcNubarn;

};

#endif
