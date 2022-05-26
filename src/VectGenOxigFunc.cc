/**
 * :file VectGenOxigFunc.cc
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
#include "VectGenOxigFunc.hh"

using namespace std;

VectGenOxigFunc::VectGenOxigFunc(){}

double VectGenOxigFunc::AngleRecCC(int num, int ix, int ex, int ch, double enu, double cos)
{
	/*
	 * calculate probability of the direction the electron is emitted after the nu_e + O charged current reaction.
	 */

	double costheta = 0;

	ifstream ifs[2][5];
	string rctn[2] = {"nue", "neb"};

	ifs[num][ix].open(Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(),ix));

	if(!ifs[num][ix]){
		cerr << "file does not exist" << endl;
		cerr << "file name is" << " " << Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(),ix) << endl;
		exit(-1);
	}

	int fileSize[2][5] = {{0}};
	int exNum[2][5] = {{3,15,8,1,6}};
	int num_ex = 0;
	double nuEne = 0.;
	std::vector<double> Enu[2][5];
	std::vector<double> totcrs[2][5];
	std::vector<double> Ee[2][5];
	int index = 0;
	double tmp_ene, tmp_crs;
	std::vector<double> nuene[2][5];
	double tmp_e0;
	double tmp;
	double tmp_cs0, tmp_cs1, tmp_cs2, tmp_cs3, tmp_cs4, tmp_cs5, tmp_cs6, tmp_sum, tmp_rec, tmp_pro, tmp_ex, tmp_num;
	std::vector<double> exEne[2][5];
	std::vector<double> rec[2][5];
	std::vector<double> pro_energy[2][5];
	if(ix==0) num_ex = 3;
	else if(ix==1) num_ex = 15;
	else if(ix==2) num_ex = 8;
	else if(ix==3) num_ex = 1;
	else if(ix==4) num_ex = 16;
	//std::cout << ix << " " << "num_ex is" << " " << num_ex << std::endl;
	while(ifs[num][ix]>>tmp>>tmp_e0){
		nuene[num][ix].push_back(tmp_e0);
		for(int j=0;j<num_ex;j++){
			ifs[num][ix]>>tmp_num>>tmp_ex>>tmp_rec>>tmp_pro>>tmp_sum>>tmp_cs0>>tmp_cs1>>tmp_cs2>>tmp_cs3>>tmp_cs4>>tmp_cs5>>tmp_cs6;
			index++;
			exEne[num][ix].push_back(tmp_ex);
			rec[num][ix].push_back(tmp_rec);
			pro_energy[num][ix].push_back(tmp_pro);
		}
	}
	fileSize[num][ix] = index;

	double me = 0.5109; // MeV

	if(ch!=8){
		double rec_energy[16] = {0.};
		rec_energy[ex] = enu - exEne[num][ix].at(ex); //energy [MeV] of recoil electron or positron

		double a = 1. + (rec_energy[ex]/25.)*(rec_energy[ex]/25.)*(rec_energy[ex]/25.)*(rec_energy[ex]/25.);
		double b = 3. + (rec_energy[ex]/25.)*(rec_energy[ex]/25.)*(rec_energy[ex]/25.)*(rec_energy[ex]/25.);

		costheta = 0.5 * (1.-(a/b)*cos);

		//return costheta;
	}

	if(ch==8){
		double recEnergy = 0.;
		if(num==0){
			recEnergy = enu - 15.4;
		}
		if(num==1){
			recEnergy = enu - 11.4;
		}
		double a = 1. + (recEnergy/25.)*(recEnergy/25.)*(recEnergy/25.)*(recEnergy/25.);
		double b = 3. + (recEnergy/25.)*(recEnergy/25.)*(recEnergy/25.)*(recEnergy/25.);

		costheta = 0.5 * (1.-(a/b)*cos);

	}
	return costheta;
}

double VectGenOxigFunc::RecEneCC(int num, int ix, int ex, int ch, double enu)
{
	/*
	 * calculate recoil e- or e+ energy 
	 */

	double recEnergy = 0;

	if(ch!=8){
		ifstream ifs[2][5];
		string rctn[2] = {"nue", "neb"};

		ifs[num][ix].open(Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(),ix));

		if(!ifs[num][ix]){
			cerr << "file does not exist" << endl;
			cerr << "file name is" << " " << Form("/disk1/disk02/usr6/nakanisi/SuperNova/time-vect300/oscillation/oxygen/crsox%s%dstate.dat",rctn[num].c_str(), ix) << endl;
			exit(-1);
		}

		int fileSize[2][5] = {{0}};
		int exNum[2][5] = {{3,15,8,1,6}};
		int num_ex = 0;
		double nuEne = 0.;
		std::vector<double> Enu[2][5];
		std::vector<double> Ee[2][5];
		int index = 0;
		double tmp_ene, tmp_crs;
		std::vector<double> nuene[2][5];
		double tmp_e0;
		double tmp;
		double tmp_cs0, tmp_cs1, tmp_cs2, tmp_cs3, tmp_cs4, tmp_cs5, tmp_cs6, tmp_sum, tmp_rec, tmp_pro, tmp_ex, tmp_num;
		std::vector<double> exEne[2][5];
		std::vector<double> rec[2][5];
		std::vector<double> pro_energy[2][5];
		if     (ix==0) num_ex = 3;
		else if(ix==1) num_ex = 15;
		else if(ix==2) num_ex = 8;
		else if(ix==3) num_ex = 1;
		else if(ix==4) num_ex = 16;
		while(ifs[num][ix]>>tmp>>tmp_e0){
			nuene[num][ix].push_back(tmp_e0);
			for(int j=0;j<num_ex;j++){
				ifs[num][ix]>>tmp_num>>tmp_ex>>tmp_rec>>tmp_pro>>tmp_sum>>tmp_cs0>>tmp_cs1>>tmp_cs2>>tmp_cs3>>tmp_cs4>>tmp_cs5>>tmp_cs6;
				index++;
				exEne[num][ix].push_back(tmp_ex);
				rec[num][ix].push_back(tmp_rec);
				pro_energy[num][ix].push_back(tmp_pro);
			}
		}
		fileSize[num][ix] = index;

		double me = 0.5109; // MeV

		double rec_energy[16] = {0.};
		rec_energy[ex] = enu - exEne[num][ix].at(ex); //energy [MeV] of recoil electron or positron

		recEnergy = rec_energy[ex];
	}

	if(ch==8){
		if(num==0){
			//if(enu > 15.4){
				recEnergy = enu - 15.4;
			//}
			//else recEnergy = 0.;
		}
		if(num==1){
			//if(enu > 11.4){
				recEnergy = enu - 11.4;
			//}
			//else recEnergy = 0.;
		}
	}

	return recEnergy;
}
