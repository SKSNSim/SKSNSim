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
#include "VectGenSetBin.hh"
#include "VectGenOxyCrosssection.hh"

using namespace std;

VectGenOxyCrosssection::VectGenOxyCrosssection()
:isOpened{false},isOpened2{false}
{}

VectGenOxyCrosssection::~VectGenOxyCrosssection()
{
  for (int num=0; num<2; num++) {
    for (int ix=0; ix<5; ix++) {
      ifs[num][ix].close();
    }
  }

  for (int num=0; num<2; num++) {
    ifs_nc[num].close();
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
				if((enu-exEne[num][ix].at(ex))>Me){
					rec_energy[ex] = enu - exEne[num][ix].at(ex);
					//std::cout << i << " " << ix << " " << num_ex << " " << j << " " << enu << " " << exEne[num][ix].at(j) << " " << rec_energy[j] << std::endl;
				}
				else{
					rec_energy[ex] = 0.;
				}
				//if(rec_energy[ex]>5.){
				if(rec_energy[ex]>eEneThr){
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
double VectGenOxyCrosssection::CsNuOxyNCNue(int num, int ix, double enu)
{

	/*
	  Total cross section of nu_e + O --> 15N + p + gamma, nu_e_bar + O --> 15O + n + gamma interaction
	*/

	double totcso = 0.;

	string rctn[2] = {"N", "O"};

	if (!isOpened2[num]) {
		std::cout << "new file open (num, ix, isOpen2): " << num << " " << isOpened2[num] << std::endl; //nakanisi
		//ifs_nc[num].open(Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/total_nc_cross_nu_%s2.dat",rctn[num].c_str()));
		ifs_nc[num].open(Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crossNC_ex%s.dat",rctn[num].c_str()));

		if(!ifs_nc[num]){
			cerr << "file does not exist" << endl;
			cerr << "file name is" << " " << Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crossNC_ex%s.dat",rctn[num].c_str()) << endl;
			exit(-1);
		}

		while(ifs_nc[num]>>tmp_nc_e0){
			nuene_nc[num].push_back(tmp_nc_e0);
			//std::cout << "neutrino energy: " << tmp_nc_e0 << std::endl;
			for(int j=0;j<num_ex_nc[num];j++){
				ifs_nc[num]>>tmp_nc_ex>>tmp_nc_cs0;
				//std::cout << num << " " << ix << " " << j << " " << index_nc[num] << " " << tmp_nc_ex << " " << tmp_nc_cs0 << std::endl;
				index_nc[num]++;
				exEne_nc[num][j].push_back(tmp_nc_ex);
				if(num==0)crs0_nc[j].push_back(tmp_nc_cs0);
				else if(num==1)crs1_nc[j].push_back(tmp_nc_cs0);
				//crs1_nc[num].push_back(tmp_cs1);
			}
		}
		isOpened2[num] = true;
	}

	//std::cout << "index_nc " << num << " " << ix << " " << index_nc[num] << std::endl;
	for(int i=0;i<20;i++){
		//if(enu>6.)std::cout << i << " " << num << " " << ix << " " << enu << " " << exEne_nc[num][ix].at(i) << std::endl;
		if(i>0 && enu<nuene_nc[num].at(i) && enu>nuene_nc[num].at(i-1) && enu>exEne_nc[num][ix].at(i)){
			if(num==0){
				totcso = ((((crs0_nc[ix].at(i)-crs0_nc[ix].at(i-1))/(nuene_nc[num].at(i)-nuene_nc[num].at(i-1)))*enu + crs0_nc[ix].at(i-1) - ((crs0_nc[ix].at(i)-crs0_nc[ix].at(i-1))/(nuene_nc[num].at(i)-nuene_nc[num].at(i-1)))*nuene_nc[num].at(i-1))) * 1.0e-42;
				//std::cout << "i " << i << " num " << num  << " ix " << ix << " enu " << enu << " exEne " << exEne_nc[num][ix].at(i) << " nuene " << nuene_nc[num].at(i) << " crs0 " << crs0_nc[ix].at(i) << " " << crs0_nc[ix].at(i-1) << " " << totcso << std::endl;
			}
			else if(num==1){
				totcso = ((((crs1_nc[ix].at(i)-crs1_nc[ix].at(i-1))/(nuene_nc[num].at(i)-nuene_nc[num].at(i-1)))*enu + crs1_nc[ix].at(i-1) - ((crs1_nc[ix].at(i)-crs1_nc[ix].at(i-1))/(nuene_nc[num].at(i)-nuene_nc[num].at(i-1)))*nuene_nc[num].at(i-1))) * 1.0e-42;
				//std::cout << "i " << i << " num " << num  << " ix " << ix << " enu " << enu << " exEne " << exEne_nc[num][ix].at(i) << " nuene " << nuene_nc[num].at(i) << " crs1 " << crs1_nc[ix].at(i) << std::endl;
			}
			break;
		}
		else{
			//std::cout << "num " << num  << " ix " << ix << " enu " << enu << " nuene " << std::endl;
			totcso = 0.;
			//break;
		}
	}

	return totcso;
}



