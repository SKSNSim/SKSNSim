/**
 * @file VectGenSnFlux.hh
 *
 * @date 2021-09-22
 * @author Y.Koshio
 */

#ifndef VECTGENSNFLUX_HH
#define VECTGENSNFLUX_HH

#include <string>
#include <vector>

/**
 * @class VectGenSnFlux
 */
class VectGenSnFlux
{

public:

  VectGenSnFlux(){}
  ~VectGenSnFlux(){}

  double VectGenSnNspeNue(double, double);
  double VectGenSnNspeNueb(double, double);
  double VectGenSnNspeNux(double, double);

protected:

  std::vector<double> tmesh;
  std::vector<std::vector<double> > enue, eneb, enux, nnue, nneb, nnux, lnue, lneb, lnux;

};

#endif
