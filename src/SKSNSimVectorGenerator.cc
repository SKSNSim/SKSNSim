/***********************************
 * SKSNSimVectorGenerator.cc
 * *********************************/
#include <functional>
#include <TRandom3.h>
#include <geotnkC.h>
#include "SKSNSimVectorGenerator.hh"
#include "SKSNSimConstant.hh"

using namespace SKSNSimPhysConst;

template<class T>
class UtilVector3 {
  public:
    T x,y,z;
    UtilVector3(T xx, T yy, T zz) : x(xx), y(yy), z(zz){}
    T operator[](size_t i)const { if(i==0) return x; if(i==1) return y; if(i==2) return z; return 9999.;}
//    UtilVector3 operator*(const T a){ return UtilVector3( a* x, a* y, a*z); }
    UtilVector3 operator+(const UtilVector3 v) const { return UtilVector3( v.x + x,  v.y+ y, v.z+z); }
    UtilVector3 operator-(const UtilVector3 v) const { return -1.0 * v + *this; }
    T Mag2() const { return x*x + y*y + z*z; }
    T Mag() const { return std::sqrt(Mag2()); }
    UtilVector3 Unit() const { return (1./Mag())* *this ; }
};
template <class T> UtilVector3<T> operator*(T a, const UtilVector3<T> &v){ return UtilVector3<T>(a*v.x, a*v.y, a*v.z);}

double SKSNSimVectorGenerator::FindMaxProb ( SKSNSimFluxModel &flux, SKSNSimCrosssectionModel &xsec){
  double maxP = 0.;
  const double ene_min = flux.GetEnergyLimitMin();
  const double ene_max = flux.GetEnergyLimitMax();
  const double cost_min = -1.;
  const double cost_max =  1.;
  const size_t nbin_ene = 1000; // TODO define appropriate binsize
  const size_t nbin_cost = 1000; // TODO define appropriate binsize
  const double diff_ene = (ene_max - ene_min)/nbin_ene;
  const double diff_cost = (cost_max - cost_min)/nbin_cost;
  for(size_t i = 0; i < nbin_ene; i++){
    const double ene =  ene_min +  diff_ene* double(i);
    for(size_t j = 0; j < nbin_cost; j++){
      const double cost =  cost_min +  diff_cost* double(i);
      const double p = flux.GetFlux(ene) * xsec.GetDiffCrosssection(ene, cost).second;
      if(maxP < p) maxP = p;
    }
  }
  return maxP;
}


SKSNSimSNEventVector SKSNSimVectorGenerator::GenerateEventIBD() {
  SKSNSimSNEventVector ev;

  double nuEne, cost, eEne;
  double nEne;

  if(fluxmodels.size() == 0) return ev;
  SKSNSimFluxModel &flux = *fluxmodels[0]; // TODO modify for user to select models
  if(xsecmodels.size() == 0) return ev;
  SKSNSimCrosssectionModel &xsec = *xsecmodels[0]; // TODO modify for user to select models

  const double maxProb = FindMaxProb(*fluxmodels[0], *xsecmodels[0]);
  auto getRandomReal = std::bind([](double min, double max, TRandom &rng){ return rng.Uniform(max - min) + min; }, std::placeholders::_1, std::placeholders::_2, randomgenerator);
  auto SQ = [](double a){ return a*a;};

  // determine neutrino and positron energy, and its direction
  while( 1 ){
    nuEne = getRandomReal( GetEnergyMin(), GetEnergyMax());

    const double nuFlux = flux.GetFlux(nuEne);

    double sigm;
    cost = getRandomReal( -1., 1.);
    auto xsecpair = xsec.GetDiffCrosssection(nuEne, cost);
    eEne = xsecpair.first;
    sigm = xsecpair.second;

    double p = nuFlux * sigm;
    double x = getRandomReal( 0., maxProb);
    if( x < p ) break;
  }

  // determine neutrino direction
  //double nuDir[3];
  double theta = getRandomReal( 0., 1.) * M_PI;
  double phi = getRandomReal( 0., 1.) * 2. * M_PI;
  UtilVector3<double> nuDir (
      std::sin(theta) * std::cos(phi),
      sin( theta ) * sin( phi ),
      cos( theta )
      );


  // Rotation matrix of neutrino direction
  double Rmat[3][3];
  Rmat[0][0] = cos( theta ) * cos( phi );
  Rmat[0][1] = -sin( phi );
  Rmat[0][2] = sin( theta ) * cos( phi );
  Rmat[1][0] = cos( theta ) * sin( phi );
  Rmat[1][1] = cos( phi );
  Rmat[1][2] = sin( theta ) * sin( phi );
  Rmat[2][0] = -sin( theta );
  Rmat[2][1] = 0.;
  Rmat[2][2] = cos( theta );

  // interaction point
  auto determinePosition = std::bind([](TRandom &rng)
  {
    const double rPositionRange = RINTK;
    const double hPositionRange = ZPINTK;
    //random inside the full tank (32.5kton)
    const double r2 = rng.Uniform(1.) * rPositionRange*rPositionRange ;
    const double r = std::sqrt( r2 );
    const double phi = rng.Uniform( 2. * M_PI);
    const double x = r * std::cos( phi );
    const double y = r * std::sin( phi );
    const double z = -hPositionRange + rng.Uniform( 2.*hPositionRange);

    return UtilVector3<double>(x,y,z);
  }, randomgenerator);
  const auto xyz = determinePosition();

  // Fill into class
  // MCVERTEX (see $SKOFL_ROOT/inc/vcvrtx.h )                                                                               
  ev.AddVertex(xyz.x, xyz.y, xyz.z,
      1, 0, 0.);


  // IBD interaction

  // Original neutrino
  const auto nuMomentum = nuEne * nuDir.Unit();
  const auto nuebarid = ev.AddTrack(
      -12, nuEne,
      nuMomentum.x, nuMomentum.y, nuMomentum.z,
      0, 0, 1, -1, 0
      );

  // Original proton
  const auto protonid = ev.AddTrack(
      2212, Mp,
      0., 0., 0.,
      0, 0, 1, -1, 0
      );

  // Positron
  double amom = sqrt(SQ( eEne ) - SQ( Me ));
  double eTheta = acos( cost );
  double ePhi = getRandomReal( -M_PI,  M_PI);

  // conversion the positron direction along the neutrino direction

  UtilVector3<double> origVec(
      sin( eTheta ) * cos( ePhi ),
      sin( eTheta ) * sin( ePhi ),
      cos( eTheta ));

  UtilVector3<double> eDir(0., 0., 0.);
  {
    double eDirBuf[3];
    eDirBuf[0] = 0.;
    eDirBuf[1] = 0.;
    eDirBuf[2] = 0.;
    for( int i = 0; i < 3; i++ ){
      for( int j = 0; j < 3; j++ ){
        eDirBuf[i] += Rmat[i][j] * origVec[j];
      }
    }
    eDir = UtilVector3<double>(eDirBuf[0], eDirBuf[1], eDirBuf[2]);
  }

  //------------------------------------------------------------------------

  // Positron
  const auto positronMomentum = amom * eDir.Unit();
  const auto positronid = ev.AddTrack(
      -11, eEne,
      positronMomentum.x, positronMomentum.y, positronMomentum.z,
      1, 1, 1, 0, 1
      );

  // Neutron
  const UtilVector3<double> neutronMomentum = nuMomentum - positronMomentum;
  const auto neutronid = ev.AddTrack(
      2112,
      std::sqrt(neutronMomentum.Mag2() + Mn*Mn),
      neutronMomentum.x, neutronMomentum.y, neutronMomentum.z,
      1, 1, 1, 0, 1
      );

  ev.SetRunnum(999999); // TODO modify
  ev.SetSubRunnum(0);  // TODO

  return ev;
}
