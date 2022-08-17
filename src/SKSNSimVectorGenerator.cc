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
    T operator*(const UtilVector3 &v) const { return v.x * x + v.y * y + v.z*z; }
    UtilVector3 operator+(const UtilVector3 v) const { return UtilVector3( v.x + x,  v.y+ y, v.z+z); }
    UtilVector3 operator-(const UtilVector3 v) const { return -1.0 * v + *this; }
    T Mag2() const { return *this * *this; }
    T Mag() const { return std::sqrt(Mag2()); }
    UtilVector3 Unit() const { return (1./Mag())* *this ; }
};
template <class T> UtilVector3<T> operator*(T a, const UtilVector3<T> &v){ return UtilVector3<T>(a*v.x, a*v.y, a*v.z);}

template<class T>
class UtilMatrix3 {
  public:
    T v[3][3];
    UtilMatrix3(){}
    ~UtilMatrix3(){}
    UtilMatrix3(
        T x00, T x01, T x02,
        T x10, T x11, T x12,
        T x20, T x21, T x22
        ){
      v[0][0] = x00; v[0][1] = x01; v[0][2] = x02;
      v[1][0] = x10; v[1][1] = x11; v[1][2] = x12;
      v[2][0] = x20; v[2][1] = x21; v[2][2] = x22;
    }
    UtilMatrix3(
        UtilVector3<T> v0, UtilVector3<T> v1, UtilVector3<T> v2
        ){
      v[0][0] = v0.x; v[0][1] = v1.x; v[0][2] = v2.x;
      v[1][0] = v0.y; v[1][1] = v1.y; v[1][2] = v2.y;
      v[2][0] = v0.z; v[2][1] = v1.z; v[2][2] = v2.z;
    }
    UtilVector3<T> operator*(const UtilVector3<T> &l) const {
      T x[3];
      for(int i = 0; i < 3; i++){
        x[i] = UtilVector3<T>(v[i][0], v[i][1], v[i][2]) * l;
      }
      return UtilVector3<T>(x[0], x[1], x[2]);
    }

    UtilMatrix3 operator*(const UtilMatrix3 &m) const {
      auto r = UtilMatrix3();
      for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
          r[i][j] = 0;
          for(int k = 0; k < 3; k++)
            r[i][j] += this->v[i][k] * m.v[k][j];
        }
      }
      return r;
    }
};

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
      const double p = flux.GetFlux(ene, 0.0, SKSNSimFluxModel::FLUXNUEB) * xsec.GetDiffCrosssection(ene, cost).first;
      if(maxP < p) maxP = p;
    }
  }
  return maxP;
}


double SKSNSimVectorGenerator::SetMaximumHitProbability(){
  // TODO fix to use appropriate combination of flux and xsec models
  if( fluxmodels.size() < 1 || xsecmodels.size() < 1 ){
    m_max_hit_probability = -1.0;
    return m_max_hit_probability;
  }
  m_max_hit_probability = FindMaxProb(*fluxmodels[0], *xsecmodels[0]);
  return m_max_hit_probability;
}

SKSNSimSNEventVector SKSNSimVectorGenerator::GenerateEventIBD() {
  SKSNSimSNEventVector ev;

  double nuEne, cost, eEne;
  double nEne;

  TRandom &rng = *randomgenerator;
  if(fluxmodels.size() == 0) return ev;
  SKSNSimFluxModel &flux = *fluxmodels[0]; // TODO modify for user to select models
  if(xsecmodels.size() == 0) return ev;
  SKSNSimCrosssectionModel &xsec = *xsecmodels[0]; // TODO modify for user to select models

  auto SQ = [](double a){ return a*a;};

  // determine neutrino and positron energy, and its direction
  while( 1 ){
    nuEne = rng.Uniform( GetEnergyMin(), GetEnergyMax());

    const double nuFlux = flux.GetFlux(nuEne, 0.0, SKSNSimFluxModel::FLUXNUEB);

    double sigm;
    cost = rng.Uniform( -1., 1.);
    auto xsecpair = xsec.GetDiffCrosssection(nuEne, cost);
    eEne = xsecpair.second;
    sigm = xsecpair.first;

    double p = nuFlux * sigm;
    double x = rng.Uniform( 0., m_max_hit_probability);
    if( x < p ) break;
  }
#ifdef DEBUG
  std::cout << "In GenerateEventIBD: eEne = " << eEne << std::endl;
#endif

  // determine neutrino direction
  //double nuDir[3];
  double theta = rng.Uniform(M_PI);
  double phi = rng.Uniform(2. * M_PI);
  UtilVector3<double> nuDir ( std::sin(theta)*std::cos(phi),
                                        sin(theta)*sin(phi),
                                                 cos(theta));


  // Rotation matrix of neutrino direction
  UtilMatrix3<double> Rmat( cos(theta)*cos(phi), -sin(phi), sin(theta)*cos(phi),
                            cos(theta)*sin(phi),  cos(phi), sin(theta)*sin(phi),
                                    -sin(theta),        0.,          cos(theta));

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
  }, *randomgenerator);
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
  double ePhi = rng.Uniform( -M_PI,  M_PI);

  // conversion the positron direction along the neutrino direction

  UtilVector3<double> origVec(
      sin( eTheta ) * cos( ePhi ),
      sin( eTheta ) * sin( ePhi ),
      cos( eTheta ));

  const UtilVector3<double> eDir = Rmat * origVec;

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
