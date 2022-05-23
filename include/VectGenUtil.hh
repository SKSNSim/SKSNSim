#ifndef VECTGENUTIL_H_INCLUDED
#define VECTGENUTIL_H_INCLUDED

#define SQ( x ) ( ( x ) * ( x ) )

const double zero_precision = 0.000001;

#include "VectGenSnRoot.hh"

//generator double random number in the range [ds, de].
double getRandomReal( const double ds, const double de , TRandom3* gen)
{
	double x = gen->Rndm();
	double p = ds + ( de - ds ) * x;
	return p;
}

#include "snevtinfo.h"
bool evtInfoTSort( const SNEvtInfo & a, const SNEvtInfo & b )
{
	return ( a.rTime < b.rTime );
}

#include <iostream>
#include <string>
#include <vector>

std::vector<std::string> split(std::string str, char del) {
  long unsigned int first = 0;
  long unsigned int last = str.find_first_of(del);
 
  std::vector<std::string> result;
  while (first < str.size()) {
    std::string subStr(str, first, last - first);
 
    result.push_back(subStr);
 
    first = last + 1;
    last = str.find_first_of(del, first);
 
    if (last == std::string::npos) {
      last = str.size();
    }
  }
  return result;
}

#endif // SETBIN_H_INCLUDED
