/**********************************
 * File: SKSNSimFlux.cc
 * Desctiption:
 ************************************/

#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>
#include "SKSNSimFlux.hh"

#ifdef RETURNEXCEPATIONS
#include <stdexcept>
#endif

SKSNSimDSNBFluxCustom::SKSNSimDSNBFluxCustom()
{
      ene_flux_v = std::make_unique<std::vector<std::pair<double,double>>>();
      lower_energy_bin_width = 0.1;
      upper_energy_bin_width = 0.1;
      return;
}

SKSNSimDSNBFluxCustom::SKSNSimDSNBFluxCustom(const std::string fname)
{
      std::cout << " DSNB flux file: "<<fname<< std::endl;
      ene_flux_v = std::make_unique<std::vector<std::pair<double,double>>>();
      this->loadFile(fname);
      return;
}

void SKSNSimDSNBFluxCustom::loadFile(const std::string fname)
{
  ene_flux_v->clear();
  std::ifstream datafile(fname.c_str());
  double ene, flux;
  for(std::string line; std::getline(datafile, line); ){
    std::stringstream ss(line);
    ss >> ene >> flux;
    ene_flux_v->push_back(std::make_pair(ene, flux));
  }
  datafile.close();
  sortByEnergy();
  return;
}

void SKSNSimDSNBFluxCustom::DumpFlux(std::ostream &out) const 
{
  for(std::vector<std::pair<double, double> >::iterator it = ene_flux_v->begin(); it != ene_flux_v->end(); it++)
    out << it->first << " MeV, " << it->second << std::endl;
}

void SKSNSimDSNBFluxCustom::sortByEnergy()
{
  if(ene_flux_v == NULL) return;
  // for std::pair, first comparing first-component then comparing second-component
  std::sort(ene_flux_v->begin(), ene_flux_v->end());
  if(ene_flux_v->size() > 1){
    lower_energy_bin_width = (ene_flux_v->begin() + 1)->first - ene_flux_v->begin()->first;
    upper_energy_bin_width = (ene_flux_v->end()-1)->first - (ene_flux_v->end()-2)->first;
  }
  return;
}

double SKSNSimDSNBFluxCustom::GetFlux(const double nu_ene_MeV, const double time_sec, const FLUXNUTYPE nutype) const 
{
  // Assumed the data field ene_flux_v is sorted as lowest energy on first 
  // Calculate flux with linear interpolation
  const static double ERROR_CODE = -9999.;
  const static double OUTOFRANGE = -9998.;
  if(ene_flux_v == NULL) return ERROR_CODE;
  std::pair<double,double> bin = std::make_pair(-1, 0);
  std::pair<double,double> nextbin;
  for(std::vector<std::pair<double,double> >::iterator it = ene_flux_v->begin(); it+1 != ene_flux_v->end(); it++){
    // IMPROVEMENT: if you implement this with binary tree search, it will be faster
    if(it->first <= nu_ene_MeV && nu_ene_MeV < (it+1)->first){
      bin = *it;
      nextbin = *(it+1);
      break;
    }
  }
  if(bin.first == -1) {
#ifdef RETURNEXCEPATIONS
    throw  std::out_of_range("energy is out of range");
#endif
    return OUTOFRANGE;
  }

  double nuFlux = (nextbin.second - bin.second) * (nu_ene_MeV - bin.first) / (nextbin.first - bin.first) + bin.second;
#ifdef DEBUG
  std::cout << "FluxCalculation: nuEne = " << nu_ene_MeV <<", flux = " << nuFlux << ", bin = (" << bin.first << ", " << bin.second << ")" << std::endl;
#endif

  return nuFlux;
}


SKSNSimSNFluxCustom::LoadFluxFile(std::string fname){
  cout <<"SN model data in SnLoading :  "<<fname << endl;

  // file open
  ifstream ifs(file_name.c_str());
  if(!ifs.is_open()){
    cerr<<"file load failed"<<endl;
    exit(-1);
  }

  // Count Energy bin from the table
  std::string line;
  int Ebin = 0;
  while(std::getline(ifs, line)){
    if(line.size() < 2) break;
    Ebin++;
  }
  Ebin--;
  ifs.close();

  // Read data
  ifs.open(file_name.c_str());
  double t0, elow, ehigh;
  vector<double> a0(Ebin), a1(Ebin), a2(Ebin), a3(Ebin), a4(Ebin), a5(Ebin), a6(Ebin), a7(Ebin), a8(Ebin);

  while(ifs>>t0){
    tmesh.push_back(t0);
    for(int j=0; j<Ebin ;j++){
      ifs>>elow>>ehigh>>a3[j]>>a4[j]>>a5[j]>>a6[j]>>a7[j]>>a8[j];

      if(a3[j] > zero_precision) a0[j] = a6[j]/a3[j]*ERG2MEV;
      else {
        a0[j] = (elow+ehigh)/2.;
        a3[j] = 0;
        a6[j] = 0;
      }
      if(a4[j] > zero_precision) a1[j] = a7[j]/a4[j]*ERG2MEV;
      else {
        a1[j] = (elow+ehigh)/2.;
        a4[j] = 0.;
        a7[j] = 0.;
      }
      if(a5[j] > zero_precision) a2[j] = a8[j]/a5[j]*ERG2MEV;
      else {
        a2[j] = (elow+ehigh)/2.;
        a5[j] = 0.;
        a8[j] = 0.;
      }
      //std::cout << t0 << " " << a0[j] << " " << a1[j] << " " << a2[j] << std::endl;// kasiwagi
      //std::cout << a3[j] << " " << a4[j] << " " << a5[j] << std::endl;// kasiwagi

    }
    enue.push_back(a0);
    eneb.push_back(a1);
    enux.push_back(a2);
    nnue.push_back(a3);
    nneb.push_back(a4);
    nnux.push_back(a5);
    lnue.push_back(a6);
    lneb.push_back(a7);
    lnux.push_back(a8);
  }

  ifs.close();

  return;
}

double SKSNSimSNFluxCustom::GetFlux(const double e, const double t, const FLUXNUTYPE type) const {

  double nspc = 0;
  double nspc0 = 0.;

  if(t < tmesh[0] || t < tmesh[tmesh.size()-1]){
    int i = 0;
    while(tmesh[i] < t) i++;
    i--;
    int j = 1;
    while(enue[i][j] < e && j<20) j++;

    //cout << "time   " << time << " " << i << " " << tmesh[i] << " " << tmesh[i+1] << endl;
    //cout << "energy " << energy << " " << j << endl;

    double nspclow, nspchigh, elow, ehigh;

    nspclow  = nnue[i][j-1];
    nspchigh = nnue[i][j];
    elow  = enue[i][j-1];
    ehigh = enue[i][j];
    if(ehigh!=0 || elow!=0){
      nspc0 = (nspchigh - nspclow) * (e - elow) / (ehigh - elow) + nspclow;
    }
    else if(ehigh==0 && elow==0){
      nspc0 = 0.;
    }
    //cout << nspc0 << endl;

    nspclow  = nnue[i+1][j-1];
    nspchigh = nnue[i+1][j];
    elow  = enue[i+1][j-1];
    ehigh = enue[i+1][j];
    if(ehigh!=0 || elow!=0){
      double nspc1 = (nspchigh - nspclow) * (e - elow) / (ehigh - elow) + nspclow;
      //cout << nspc1 << endl;

      nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;
    }
    else if(ehigh==0 && elow==0){
      double nspc1 = 0.;
      nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;
    }
    //cout << nspc << endl;
  }

  return nspc;

}
