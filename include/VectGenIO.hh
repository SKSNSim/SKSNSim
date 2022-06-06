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
}

class VectGenIO : VectGenGenerator
{

public:
  VectGenIO(){}
  VectGenIO(std::string, int, double, int, std::string, uint);
  VectGenIO(std::string, uint, std::string, double emin = -1.0, double emax = -1.0); // For DSBN vector generator
  ~VectGenIO(){}

  void DoProcess();
  void DoProcess(int);

private:

  std::string FileText;

};

#endif
