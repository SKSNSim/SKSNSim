#ifndef __SKSNSIMTOOLS_H_INCLUDED__
#define __SKSNSIMTOOLS_H_INCLUDED__

#include <TString.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "SKSNSimEnum.hh"

namespace SKSNSimTools{
  auto DumpDebugMessage = [] (TString str) {
#ifdef DEBUG
    std::cout << "[IZU:DEBUG] " << str << std::endl;
#endif
  };
}

namespace SKSNSimLiveTime {
  
  std::vector<std::tuple<int /* runnum */,double /* livetime_day */>> LoadLiveTime(SKSNSIMENUM::SKPERIOD);
  std::vector<std::tuple<int /* runnum */,double /* livetime_day */>> LoadLiveTime(int /* run-begin */,int /* run-end */);
  std::vector<std::tuple<int /* runnum */,double /* livetime_day */>> LoadLiveTime(std::string);

}
#endif
