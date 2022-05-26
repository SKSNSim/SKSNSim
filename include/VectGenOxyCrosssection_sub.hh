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

	double CsNuOxy43CCSub(int, int, int, double);
};

#endif
