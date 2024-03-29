/**
 * @class Crosssection
 */

#ifndef VECTGENOXYCROSSSECTION_HH
#define VECTGENOXYCROSSSECTION_HH
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"


class VectGenOxyCrosssection
{
public:
	VectGenOxyCrosssection();
	~VectGenOxyCrosssection();

	double CsNuOxy43CC(int, int, int, int, double);

	double CsNuOxyNCNue(int, int, double);
	//double CsNuOxyNCNue(int, double);

private:
	bool isOpened[2][5];
	bool isOpened2[2];
	ifstream ifs[2][5];
	ifstream ifs_nc[2];
	std::vector<double> exEne[2][5];
	std::vector<double> rec[2][5];
	std::vector<double> pro_energy[2][5];
	std::vector<double> crs0[2][5][16];
	std::vector<double> crs1[2][5][16];
	std::vector<double> crs2[2][5][16];
	std::vector<double> crs3[2][5][16];
	std::vector<double> crs4[2][5][16];
	std::vector<double> crs5[2][5][16];
	std::vector<double> crs6[2][5][16];

	std::vector<double> exEne_nc[2][14];
	std::vector<double> crs0_nc[14];
	std::vector<double> crs1_nc[6];
    //std::vector<double> crs_nc[2];

	int fileSize[2][5] = {{0}};
	int exNum[2][5] = {{3,15,8,1,6}};
	double nuEne = 0.;
	std::vector<double> Enu[2][5];
	std::vector<double> totcrs[2][5];
	std::vector<double> Ee[2][5];
	int index = 0;
	int index_nc[2] = {0};
	double tmp_ene, tmp_crs;
	std::vector<double> nuene[2][5];
	std::vector<double> nuene_nc[2];
	double tmp_e0, tmp_nc_e0;
	double tmp;
	double tmp_cs0, tmp_cs1, tmp_cs2, tmp_cs3, tmp_cs4, tmp_cs5, tmp_cs6, tmp_sum, tmp_rec, tmp_pro, tmp_ex, tmp_num;
	double tmp_nc_ex, tmp_nc_cs0;
	const int num_ex[5] = {3, 15, 8, 1, 16};
	//const int num_ex_nc[2] = {14, 6};
	const int num_ex_nc[2] = {8, 4};

};

#endif
