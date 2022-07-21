/***********************************
 * SKSNSimVectorGenerator.cc
 * *********************************/
#include <functional>
#include <TRandom3.h>
#include <geotnkC.h>
#include "SKSNSimVectorGenerator.hh"
#include "SKSNSimConstant.hh"

using namespace SKSNSimPhysConst;

double FindMaxProb ( SKSNSimFluxModel &flux, SKSNSimCrosssectionModel &xsec){
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
  if(xsecmodels.size() == 0) return ev;

  const double maxProb = FindMaxProb(*fluxmodels[0], *xsecmodels[0]);
  auto getRandomReal = std::bind([](double min, double max, TRandom &rng){ return rng.Uniform(max - min) + min; }, std::placeholders::_1, std::placeholders::_2, randomgenerator);
  auto SQ = [](double a){ return a*a;};
  // determine neutrino and positron energy, and its direction                                                              
  while( 1 ){
    nuEne = getRandomReal( GetEnergyMin(), GetEnergyMax());

    const double nuFlux = fluxmodels[0]->GetFlux(nuEne);

    double sigm;
    cost = getRandomReal( -1., 1.);
    auto xsecpair = xsecmodels[0]->GetDiffCrosssection(nuEne, cost);
    eEne = xsecpair.first;
    sigm = xsecpair.second;

    double p = nuFlux * sigm;
    double x = getRandomReal( 0., maxProb);
    if( x < p ) break;
  }

  // determine neutrino direction
  double nuDir[3];
  double theta = getRandomReal( 0., 1.) * M_PI;
  double phi = getRandomReal( 0., 1.) * 2. * M_PI;
  nuDir[0] = sin( theta ) * cos( phi );
  nuDir[1] = sin( theta ) * sin( phi );
  nuDir[2] = cos( theta );

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

    return std::make_tuple(x,y,z);
  }, randomgenerator);
  const auto xyz = determinePosition();
  const double ver_x = std::get<0>(xyz);
  const double ver_y = std::get<1>(xyz);
  const double ver_z = std::get<2>(xyz);

  // Fill into class
  // MCVERTEX (see $SKOFL_ROOT/inc/vcvrtx.h )                                                                               
  ev.AddVertex(ver_x, ver_y, ver_z,
      1, 0, 0.);


  // IBD interaction

  // Original neutrino
  const auto nuebarid = ev.AddTrack(
      -12, nuEne,
      nuEne * nuDir[0], nuEne * nuDir[1], nuEne * nuDir[2],
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
  double origVec[3], eDir[3];

  origVec[0] = sin( eTheta ) * cos( ePhi );
  origVec[1] = sin( eTheta ) * sin( ePhi );
  origVec[2] = cos( eTheta );

  eDir[0] = 0.;
  eDir[1] = 0.;
  eDir[2] = 0.;

  for( int i = 0; i < 3; i++ ){
    for( int j = 0; j < 3; j++ ){
      eDir[i] += Rmat[i][j] * origVec[j];
    }
  }
  //------------------------------------------------------------------------

  // Positron
  const auto positronid = ev.AddTrack(
      -11, eEne,
      amom * eDir[0], amom * eDir[1], amom * eDir[2],
      1, 1, 1, 0, 1
      );

  // Neutron
  const double neutron_momentum[3] = {
    ev.GetTrackMomentumX(nuebarid) - ev.GetTrackMomentumX(positronid),
    ev.GetTrackMomentumY(nuebarid) - ev.GetTrackMomentumY(positronid),
    ev.GetTrackMomentumZ(nuebarid) - ev.GetTrackMomentumZ(positronid)
  };
  ev.AddTrack(
      2112,
      std::sqrt(neutron_momentum[0]*neutron_momentum[0]+neutron_momentum[1]*neutron_momentum[1]+neutron_momentum[2]*neutron_momentum[2] + Mn*Mn),
      neutron_momentum[0], neutron_momentum[1], neutron_momentum[2],
      1, 1, 1, 0, 1
      );

  ev.SetRunnum(999999); // TODO modify
  ev.SetSubRunnum(0);  // TODO

  return ev;
}
