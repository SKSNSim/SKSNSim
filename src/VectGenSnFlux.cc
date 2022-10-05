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
  double nspc0 = 0.;

  if(time < tmesh[0] || time < tmesh[tmesh.size()-1]){ // TODO is it fine ? time > tmesh[0]?
    int i = 0;
    while(tmesh[i] < time) i++;
    i--;
    int j = 1;
    while(enue[i][j] < energy && j<20) j++;

    //cout << "time   " << time << " " << i << " " << tmesh[i] << " " << tmesh[i+1] << endl;
    //cout << "energy " << energy << " " << j << endl;

    double nspclow, nspchigh, elow, ehigh;

    nspclow  = nnue[i][j-1];
    nspchigh = nnue[i][j];
    elow  = enue[i][j-1];
    ehigh = enue[i][j];
    if(ehigh!=0 || elow!=0){
    	nspc0 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
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
	double nspc1 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
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

double VectGenSnFlux::VectGenSnNspeNueb(double time, double energy)
{
  double nspc = 0.;
  double nspc0 = 0.;
	
  if(tmesh[0] < time && time < tmesh[tmesh.size()-1]){
    int i = 0;
    while(tmesh[i] < time) i++;
    i--;

    int j = 1;
    while(eneb[i][j] < energy && j<20) j++;

    //if(time>=0.0225 && time<0.023)std::cout << "vector size " << eneb.size() << " " << tmesh.size() << " energy " << energy << " " << "eneb[i][j] " << eneb[i][j] << " " << "i " << i << " " << "j " << j << std::endl;
	
    double nspclow, nspchigh, elow, ehigh;
		
    nspclow  = nneb[i][j-1];
    nspchigh = nneb[i][j];

    elow  = eneb[i][j-1];
    ehigh = eneb[i][j];

    if(ehigh!=0 || elow!=0){
	nspc0 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
    }
    else if(ehigh==0 && elow==0){
	nspc0 = 0.;
    }
	
    nspclow  = nneb[i+1][j-1];
    nspchigh = nneb[i+1][j];

    elow  = eneb[i+1][j-1];
    ehigh = eneb[i+1][j];

    if(ehigh!=0 || elow!=0){
	double nspc1 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
    	nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;
    	//if(time>=0.0055 && time<0.006)std::cout << "non zero " << nspc1 << "nspchigh " << nspchigh << " " << "nspclow " << nspclow << " " << "elow " << elow << " " << "ehigh " << ehigh << " " << "nspclow " << nspclow << " " << "energy " << energy << std::endl; //nakanisi
    }
    else if(ehigh==0 && elow==0){
	double nspc1 = 0.;
    	nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;
    	//if(time>=0.0055 && time<0.006)std::cout << "zero " << nspc1 << " " << "nspchigh " << nspchigh << " " << "nspclow " << nspclow << " " << "elow " << elow << " " << "ehigh " << ehigh << " " << "nspclow " << nspclow << " " << "energy " << energy << std::endl; //nakanisi
    }
    //if(time>=0.0055 && time<0.006)std::cout << "nspchigh " << nspchigh << " " << "nspclow " << nspclow << " " << "elow " << elow << " " << "ehigh " << ehigh << " " << "nspclow " << nspclow << " " << "energy " << energy << std::endl; //nakanisi
	
    //nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;
    //if(time>=0.0055 && time<0.006)std::cout << time << " " << energy << " " << "nspc1 " << nspc1 << " " << "ncpc0 " << nspc0 << " " << "tmesh[i] " << tmesh[i] << " " << "tmesh[i+1] " << tmesh[i+1] << " " << "nspc " << nspc << std::endl; //nakanisi
  }
	
  return nspc;
}

double VectGenSnFlux::VectGenSnNspeNux(double time, double energy)
{
  double nspc = 0;
  double nspc0 = 0.;

  if(tmesh[0] < time && time < tmesh[tmesh.size()-1]){
    int i = 0;
    while(tmesh[i] < time) i++;
    i--;
    int j = 1;
    while(enux[i][j] < energy && j<20) j++;

    //cout << "time   " << time << " " << i << " " << tmesh[i] << " " << tmesh[i+1] << endl;
    //if(time>=0.0225 && time<0.023)std::cout << "energy " << energy << " " << "enux[i][j] " << enux[i][j] << " " << "i " << i << " " << "j " << j << std::endl;

    double nspclow, nspchigh, elow, ehigh;

    nspclow  = nnux[i][j-1];
    nspchigh = nnux[i][j];
    elow  = enux[i][j-1];
    ehigh = enux[i][j];
    if(ehigh!=0 || elow!=0){
	nspc0 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
    }
    else if(ehigh==0 && elow==0){
	nspc0 = 0.;
    }
    //cout << nspc0 << endl;

    nspclow  = nnux[i+1][j-1];
    nspchigh = nnux[i+1][j];
    elow  = enux[i+1][j-1];
    ehigh = enux[i+1][j];
    if(ehigh!=0 || elow!=0){
    	double nspc1 = (nspchigh - nspclow) * (energy - elow) / (ehigh - elow) + nspclow;
    	nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;
    	//if(time>=0.0225 && time<0.023)std::cout << "non zero " << "i " << i << " " << "j " << j << " " << nspc1 << " " << "nspchigh " << nspchigh << " " << "nspclow " << nspclow << " " << "elow " << elow << " " << "ehigh " << ehigh << " " << "nspclow " << nspclow << " " << "energy " << energy << std::endl; //nakanisi
    }
    else if(ehigh==0 && elow==0){
	double nspc1=0.;
    	nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;
    	//if(time>=0.0225 && time<0.023)std::cout << "zero " << nspc1 << " " << "nspchigh " << nspchigh << " " << "nspclow " << nspclow << " " << "elow " << elow << " " << "ehigh " << ehigh << " " << "nspclow " << nspclow << " " << "energy " << energy << std::endl; //nakanisi
    }
    //cout << nspc1 << endl;

    //nspc = (nspc1 - nspc0) * (time - tmesh[i]) / (tmesh[i+1] - tmesh[i]) + nspc0;

    //cout << nspc << endl;
  }

  return nspc;
}
