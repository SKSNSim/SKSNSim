/**
 * @file VectGenSnFlux.cc
 *
 * @date 2017-10-08
 * @author Y.Koshio
 */
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "VectGenSnFlux.hh"

double VectGenSnFlux::VectGenSnNspeNue(double time, double energy)
{

  double nspc = 0;

  if(time < tmesh[0] || time < tmesh[tmesh.size()-1]){
    int i = 0;
    while(tmesh[i] < time) i++;
    i--;
    int j = 1;
    while(enue[i][j] < energy) j++;

    //cout << "time   " << time << " " << i << " " << tmesh[i] << " " << tmesh[i+1] << endl;
    //cout << "energy " << energy << " " << j << endl;

    double nspclow, nspchigh, elow, ehigh;

    nspclow  = nnue[i][j-1];
    nspchigh = nnue[i][j];
    elow  = enue[i][j-1];
    ehigh = enue[i][j];
    double nspc0 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
    //cout << nspc0 << endl;

    nspclow  = nnue[i+1][j-1];
    nspchigh = nnue[i+1][j];
    elow  = enue[i+1][j-1];
    ehigh = enue[i+1][j];
    double nspc1 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
    //cout << nspc1 << endl;

    nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;

    //cout << nspc << endl;
  }

  return nspc;
}

double VectGenSnFlux::VectGenSnNspeNueb(double time, double energy)
{
  double nspc = 0.;
	
  if(tmesh[0] < time && time < tmesh[tmesh.size()-1]){
    int i = 0;
    while(tmesh[i] < time) i++;
    i--;

    int j = 1;
    while(eneb[i][j] < energy) j++;
	
    double nspclow, nspchigh, elow, ehigh;
		
    nspclow  = nneb[i][j-1];
    nspchigh = nneb[i][j];

    elow  = eneb[i][j-1];
    ehigh = eneb[i][j];

    double nspc0 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
	
    nspclow  = nneb[i+1][j-1];
    nspchigh = nneb[i+1][j];

    elow  = eneb[i+1][j-1];
    ehigh = eneb[i+1][j];

    double nspc1 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
	
    nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;
  }
	
  return nspc;
}

double VectGenSnFlux::VectGenSnNspeNux(double time, double energy)
{
  double nspc = 0;

  if(tmesh[0] < time && time < tmesh[tmesh.size()-1]){
    int i = 0;
    while(tmesh[i] < time) i++;
    i--;
    int j = 1;
    while(enux[i][j] < energy) j++;

    //cout << "time   " << time << " " << i << " " << tmesh[i] << " " << tmesh[i+1] << endl;
    //cout << "energy " << energy << " " << j << endl;

    double nspclow, nspchigh, elow, ehigh;

    nspclow  = nnux[i][j-1];
    nspchigh = nnux[i][j];
    elow  = enux[i][j-1];
    ehigh = enux[i][j];
    double nspc0 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
    //cout << nspc0 << endl;

    nspclow  = nnux[i+1][j-1];
    nspchigh = nnux[i+1][j];
    elow  = enux[i+1][j-1];
    ehigh = enux[i+1][j];
    double nspc1 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
    //cout << nspc1 << endl;

    nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;

    //cout << nspc << endl;
  }

  return nspc;
}
