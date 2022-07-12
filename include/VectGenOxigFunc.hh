/**
 * @class Oxigfunc
 */

#ifndef VECTGENOXIGFUNC_HH
#define VECTGENOXIGFUNC_HH
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"


class VectGenOxigFunc
{
public:
	VectGenOxigFunc();

	double AngleRecCC(int, int, int, int, double, double);
	double RecEneCC(int, int, int, int, double);
};

#endif
