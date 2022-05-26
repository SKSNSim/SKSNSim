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

	double CsNuOxy43CC(int, int, int, int, double);
};

#endif
