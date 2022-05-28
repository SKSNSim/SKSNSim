/**
 * :file VectGenOxyCrosssection.cc
 *
 */
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include <math.h>

#include "VectGenSnConst.hh"
#include "VectGenOxyCrosssection.hh"

using namespace std;

VectGenOxyCrosssection::VectGenOxyCrosssection()
:isOpened{false}
{}

VectGenOxyCrosssection::~VectGenOxyCrosssection()
{
  for (int num=0; num<2; num++) {
    for (int ix=0; ix<5; ix++) {
      ifs[num][ix].close();
    }
  }
}
double VectGenOxyCrosssection::CsNuOxy43CC(int num, int ix, int ex, int ch, double enu)
{
	/*
	 * Total cross section of nu_e + O --> e^- + X, nu_e_bar + O --> e^+ + X interaction
	 */

	double totcso = 0;

	string rctn[2] = {"nue", "neb"};

	if (!isOpened[num][ix]) {
		std::cout << "new file open (num, ix, isOpen): " << num << " " << ix << " " << isOpened[num][ix] << std::endl; //nakanisi
		ifs[num][ix].open(Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(),ix));

		if(!ifs[num][ix]){
			cerr << "file does not exist" << endl;
			cerr << "file name is" << " " << Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(),ix) << endl;
			exit(-1);
		}

		//std::cout << ix << " " << "num_ex is" << " " << num_ex << std::endl;
		while(ifs[num][ix]>>tmp>>tmp_e0){
			nuene[num][ix].push_back(tmp_e0);
			for(int j=0;j<num_ex[ix];j++){
				ifs[num][ix]>>tmp_num>>tmp_ex>>tmp_rec>>tmp_pro>>tmp_sum>>tmp_cs0>>tmp_cs1>>tmp_cs2>>tmp_cs3>>tmp_cs4>>tmp_cs5>>tmp_cs6;
				index++;
				exEne[num][ix].push_back(tmp_ex);
				rec[num][ix].push_back(tmp_rec);
				pro_energy[num][ix].push_back(tmp_pro);
				crs0[num][ix][j].push_back(tmp_cs0);
				crs1[num][ix][j].push_back(tmp_cs1);
				crs2[num][ix][j].push_back(tmp_cs2);
				crs3[num][ix][j].push_back(tmp_cs3);
				crs4[num][ix][j].push_back(tmp_cs4);
				crs5[num][ix][j].push_back(tmp_cs5);
				crs6[num][ix][j].push_back(tmp_cs6);
			}
		}
		fileSize[num][ix] = index;
		//std::cout << fileSize[num][ix] << std::endl;
		//std::cout << num << " " << ix << " " << index << std::endl;
		//std::cout << exEne[num][ix].at(0) << " " << exEne[num][ix].at(1) << " " << exEne[num][ix].at(2) << std::endl;
		//
		isOpened[num][ix] = true;
	}

	double me = 0.5109; // MeV

	double rec_energy[16] = {0.};
	//double totcsnueo[16][7] = {{0.}};
	double nueEne[21] = {0., 12., 15., 18., 21., 24., 27., 30., 33., 36., 39., 42., 45., 50., 55., 60., 65., 70., 80., 90., 100.};
	for(int i=0;i<20;i++){
		//for(int j=0;j<num_ex;j++){
		double totcsnueo[16][7] = {{0.}};
			if(i > 0 && enu < nuene[num][ix].at(i) && enu > nuene[num][ix].at(i-1)){
				if((enu-exEne[num][ix].at(ex))>0.){
					rec_energy[ex] = enu - exEne[num][ix].at(ex);
					//std::cout << i << " " << ix << " " << num_ex << " " << j << " " << enu << " " << exEne[num][ix].at(j) << " " << rec_energy[j] << std::endl;
				}
				else{
					rec_energy[ex] = 0.;
				}
				if(rec_energy[ex]>5.){
						//int ie = int(enu) - int(nuene[num][ix].at(i)-nuene[num][ix].at(i-1));
						if(ch==0)totcsnueo[ex][ch] = (((crs0[num][ix][ex].at(i)-crs0[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*enu + crs0[num][ix][ex].at(i-1) - ((crs0[num][ix][ex].at(i)-crs0[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*nuene[num][ix].at(i-1)) * 1.0e-26;
						if(ch==1)totcsnueo[ex][ch] = (((crs1[num][ix][ex].at(i)-crs1[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*enu + crs1[num][ix][ex].at(i-1) - ((crs1[num][ix][ex].at(i)-crs1[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*nuene[num][ix].at(i-1)) * 1.0e-26;
						if(ch==2)totcsnueo[ex][ch] = (((crs2[num][ix][ex].at(i)-crs2[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*enu + crs2[num][ix][ex].at(i-1) - ((crs2[num][ix][ex].at(i)-crs2[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*nuene[num][ix].at(i-1)) * 1.0e-26;
						if(ch==3)totcsnueo[ex][ch] = (((crs3[num][ix][ex].at(i)-crs3[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*enu + crs3[num][ix][ex].at(i-1) - ((crs3[num][ix][ex].at(i)-crs3[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*nuene[num][ix].at(i-1)) * 1.0e-26;
						if(ch==4)totcsnueo[ex][ch] = (((crs4[num][ix][ex].at(i)-crs4[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*enu + crs4[num][ix][ex].at(i-1) - ((crs4[num][ix][ex].at(i)-crs4[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*nuene[num][ix].at(i-1)) * 1.0e-26;
						if(ch==5)totcsnueo[ex][ch] = (((crs5[num][ix][ex].at(i)-crs5[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*enu + crs5[num][ix][ex].at(i-1) - ((crs5[num][ix][ex].at(i)-crs5[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*nuene[num][ix].at(i-1)) * 1.0e-26;
						if(ch==6)totcsnueo[ex][ch] = (((crs6[num][ix][ex].at(i)-crs6[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*enu + crs6[num][ix][ex].at(i-1) - ((crs6[num][ix][ex].at(i)-crs6[num][ix][ex].at(i-1))/(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)))*nuene[num][ix].at(i-1)) * 1.0e-26;
						//if(ix==0 && ex==0 && ch==0){
						//	std::cout << ix << " " << ex << " " << ch << " " << i << " " << enu << " " << nuene[num][ix].at(i) << " " << ie << " " << int(nuene[num][ix].at(i)-nuene[num][ix].at(i-1)) << " " << nuene[num][ix].at(i) << " " << nuene[num][ix].at(i-1) << " " << crs0[num][ex].at(i-1) << " " << enu-double(ie+(nuene[num][ix].at(i)-nuene[num][ix].at(i-1))) << " " << totcsnueo[ex][ch] << std::endl;
						//}
						totcso = totcsnueo[ex][ch];
						if(totcsnueo[ex][ch] > 0.){
							break;
						}
						//totcso = totcsnueo[ex][ch];
						break;
				}
				else if(rec_energy[ex]==0){
					totcsnueo[ex][ch] = 0.;
					totcso = totcsnueo[ex][ch];
					break;
				}
				break;
			}
			else{
				totcsnueo[ex][ch] = 0.;
				totcso = totcsnueo[ex][ch];
			}
		//}
	}

	//totcso = totcsnueo[ex][ch];
	//if(ix==0 && ex==0 && ch==0)std::cout << num << " " << ix << " " << ex << " " << ch << " " << enu << " " << totcso << std::endl;
	return totcso;
}
