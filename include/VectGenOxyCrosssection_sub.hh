/**
 * @class Crosssection
 */

#ifndef VECTGENOXYCROSSSECTION_SUB_HH
#define VECTGENOXYCROSSSECTION_SUB_HH
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"


class VectGenOxyCrosssectionSub
{
public:
	VectGenOxyCrosssectionSub();
	~VectGenOxyCrosssectionSub();
		
	double CsNuOxy43CCSub(int, int, int, double);

private:
	ifstream ifs[2][5];
	bool isOpened[2][5];

	int fileSize[2][5] = {{0}};

	std::vector<double> Enu[2][5];
	std::vector<double> totcrs[2][5];
	std::vector<double> crs[2][5][32];

};

#endif
