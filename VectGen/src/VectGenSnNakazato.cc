/**
 * @file sn_nakazato.cc
 *
 * @date 2017-10-08
 * @author Y.Koshio
 */
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "VectGenSnNakazato.hh"

using namespace std;

//static const int Ebin=25; // change
static const double erg2mev=1.0e0/1.60217733e-6;

VectGenSnNakazato::VectGenSnNakazato(string file_name)
{

  const double zero_precision = 0.000001;

  cout <<"SN model data in SnLoading :  "<<file_name << endl;

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

      if(a3[j] > zero_precision) a0[j] = a6[j]/a3[j]*erg2mev;
      else {
	a0[j] = (elow+ehigh)/2.;
	a3[j] = 0;
	a6[j] = 0;
      }
      if(a4[j] > zero_precision) a1[j] = a7[j]/a4[j]*erg2mev;
      else {
	a1[j] = (elow+ehigh)/2.;
	a4[j] = 0.;
	a7[j] = 0.;
      }
      if(a5[j] > zero_precision) a2[j] = a8[j]/a5[j]*erg2mev;
      else {
	a2[j] = (elow+ehigh)/2.;
	a5[j] = 0.;
	a8[j] = 0.;
      }

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
