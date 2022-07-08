/**
 * @file VectGenIO.hh
 *
 * @date 2021-09-22
 * @author Y.Koshio
 */

#ifndef VECTGENIO_HH
#define VECTGENIO_HH

#include "VectGenSnRoot.hh"
#include "VectGenGenerator.hh"

#include "VectGenSnNakazato.hh"

extern "C" {
	void sn_sundir_( int *, int *, float *, float *, float *);
  void read_timevent_( int *, int *);
}

enum SKGEOM
{
    gSK_IV = 4,
    gSK_V  = 5,
    gSK_VI = 6,
};



class VectGenIO : VectGenGenerator
{

public:
  VectGenIO(){}
  VectGenIO(std::string, int, double, int, std::string, uint);
  VectGenIO(uint); // For DSBN vector generator
  //VectGenIO(std::string, uint, std::string); // For DSBN vector generator
  ~VectGenIO(){}

  void DoProcess();
  void DoProcess(int);

  void SetFluxFile(std::string);
  void OpenOutputFile(std::string);
  void CloseOutputFile();

  void SetRefRunNumber(int runNum) { fRefRunNum = runNum; }
  void SetNuEnergyRange(double ene_min, double ene_max) { nuEne_min = ene_min; nuEne_max = ene_max; }
  void SetUseTimeEvent(bool b) { bIsUseTimeEvent = b; }
  void SetUseFlatFlux(bool b) { bUseFlatFlux = b; }


private:

  std::string FileText;

};

#endif
