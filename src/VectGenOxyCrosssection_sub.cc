/**
 * :file totaccrs43_sub.cc
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
#include "VectGenOxyCrosssection_sub.hh"

using namespace std;

VectGenOxyCrosssectionSub::VectGenOxyCrosssectionSub(){}

double VectGenOxyCrosssectionSub::CsNuOxy43CCSub(int num, int ix, int ch, double enu)
{
	/*
	 * Total cross section of nu_e + O --> e^- + X, nu_e_bar + O --> e^+ + X interaction
	 */

	double totcso = 0.;

	ifstream ifs[2][5];
	string rctn[2] = {"nue", "neb"};

	ifs[num][ix].open(Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/sub3-1/dat/sub%s%dstate.dat",rctn[num].c_str(),ix));

	if(!ifs[num][ix]){
		cerr << "file does not exist" << endl;
		cerr << "file name is" << " " << Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/sub3-1/dat/sub%s%dstate.dat",rctn[num].c_str(),ix) << endl;
		exit(-1);
	}

	int fileSize[2][5] = {{0}};
	int exNum[2][5] = {{3,15,8,1,16}};
	int num_ex = 0;
	double nuEne = 0;
	std::vector<double> Enu[2][5];
	std::vector<double> totcrs[2][5];
	std::vector<double> crs[32];
	int index = 0;
	double tmp_e0;
	double tmp_cs0, tmp_cs1, tmp_cs2, tmp_cs3, tmp_cs4, tmp_cs5, tmp_cs6, tmp_cs7, tmp_cs8, tmp_cs9, tmp_cs10, tmp_cs11, tmp_cs12, tmp_cs13, tmp_cs14, tmp_cs15, tmp_cs16, tmp_cs17, tmp_cs18, tmp_cs19, tmp_cs20, tmp_cs21, tmp_cs22, tmp_cs23, tmp_cs24, tmp_cs25, tmp_cs26, tmp_cs27, tmp_cs28, tmp_cs29, tmp_cs30, tmp_cs31, tmp_sum;
	while(!ifs[num][ix].eof()){
		ifs[num][ix]>>tmp_e0>>tmp_cs0>>tmp_cs1>>tmp_cs2>>tmp_cs3>>tmp_cs4>>tmp_cs5>>tmp_cs6>>tmp_cs7>>tmp_cs8>>tmp_cs9>>tmp_cs10>>tmp_cs11>>tmp_cs12>>tmp_cs13>>tmp_cs14>>tmp_cs15>>tmp_cs16>>tmp_cs17>>tmp_cs18>>tmp_cs19>>tmp_cs20>>tmp_cs21>>tmp_cs22>>tmp_cs23>>tmp_cs24>>tmp_cs25>>tmp_cs26>>tmp_cs27>>tmp_cs28>>tmp_cs29>>tmp_cs30>>tmp_cs31;
		index++;
		Enu[num][ix].push_back(tmp_e0);
		crs[0].push_back(tmp_cs0);
		crs[1].push_back(tmp_cs1);
		crs[2].push_back(tmp_cs2);
		crs[3].push_back(tmp_cs3);
		crs[4].push_back(tmp_cs4);
		crs[5].push_back(tmp_cs5);
		crs[6].push_back(tmp_cs6);
		crs[7].push_back(tmp_cs7);
		crs[8].push_back(tmp_cs8);
		crs[9].push_back(tmp_cs9);
		crs[10].push_back(tmp_cs10);
		crs[11].push_back(tmp_cs11);
		crs[12].push_back(tmp_cs12);
		crs[13].push_back(tmp_cs13);
		crs[14].push_back(tmp_cs14);
		crs[15].push_back(tmp_cs15);
		crs[16].push_back(tmp_cs16);
		crs[17].push_back(tmp_cs17);
		crs[18].push_back(tmp_cs18);
		crs[19].push_back(tmp_cs19);
		crs[20].push_back(tmp_cs20);
		crs[21].push_back(tmp_cs21);
		crs[22].push_back(tmp_cs22);
		crs[23].push_back(tmp_cs23);
		crs[24].push_back(tmp_cs24);
		crs[25].push_back(tmp_cs25);
		crs[26].push_back(tmp_cs26);
		crs[27].push_back(tmp_cs27);
		crs[28].push_back(tmp_cs28);
		crs[29].push_back(tmp_cs29);
		crs[30].push_back(tmp_cs30);
		crs[31].push_back(tmp_cs31);
		//if(num==0 && ix==0)std::cout << tmp_e0 << " " << tmp_cs0 << " " << tmp_cs1 << " " << tmp_cs2 << " " << tmp_cs3 << " " << tmp_cs4 << " " << tmp_cs5 << " " << tmp_cs6 << " " << tmp_cs7 << " " << tmp_cs8 << " " << tmp_cs9 << " " << tmp_cs10 << " " << tmp_cs11 << " " << tmp_cs12 << " " << tmp_cs13 << " " << tmp_cs14 << " " << tmp_cs15 << " " << tmp_cs16 << " " << tmp_cs17 << " " << tmp_cs18 << " " << tmp_cs19 << " " << tmp_cs20 << " " << tmp_cs21 << " " << tmp_cs22 << " " << tmp_cs23 << " " << tmp_cs24 << " " << tmp_cs25 << " " << tmp_cs26 << " " << tmp_cs27 << " " << tmp_cs28 << " " << tmp_cs29 << " " << tmp_cs30 << " " << tmp_cs31 << std::endl;
	}
	fileSize[num][ix] = index-1;

	double me = 0.5109; // MeV
	double rec_energy = 0.;
	for(int i=0;i<fileSize[num][ix];i++){
		double totcsnueo[32] = {0.};
		if(i>0 && enu < Enu[num][ix].at(i) && enu > Enu[num][ix].at(i-1)){
			if(num==0){
				rec_energy = enu - 15.4;
			}
			if(num==1){
				rec_energy = enu - 11.4;
			}
			//else{
			//	rec_energy = 0.;
			//}
			if(rec_energy>0.){
				totcsnueo[ch] = (((crs[ch].at(i)-crs[ch].at(i-1))/(Enu[num][ix].at(i)-Enu[num][ix].at(i-1)))*enu + crs[ch].at(i-1) - ((crs[ch].at(i)-crs[ch].at(i-1))/(Enu[num][ix].at(i)-Enu[num][ix].at(i-1)))*Enu[num][ix].at(i-1));
				//std::cout << fileSize[num][ix] << " " << num << " " << ix << " " << ch << " " << i << " " << enu << " " << Enu[num][ix].at(i) << " " << Enu[num][ix].at(i-1) << " " << totcsnueo[ch] << std::endl;
				if(totcsnueo[ch]>0.){
					totcso = totcsnueo[ch];
					break;
				}
				if(totcsnueo[ch]<0.){
					totcso = 0.;
					break;
				}
			}
			else if(rec_energy<=0.){
				totcso = 0.;
				break;
			}
		}
		else{
			totcso = 0.;
		}
	}

	return totcso;
}
	

