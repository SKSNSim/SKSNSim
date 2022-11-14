#ifndef __SKSNSIMTOOLS_H_INCLUDED__
#define __SKSNSIMTOOLS_H_INCLUDED__

#include <TString.h>
#include <TRandom.h>
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

  SKSNSIMENUM::SKPERIOD FindSKPeriod(int /* run */);

  int elapseday(int /* year */, int /* month */, int /* day */);
  int elapseday(int /* run */);
}

namespace SKSNSimLiveTime {
  typedef std::tuple<int /* runnum */, double /*livetime_day */> RECORDLIVETIME;
  
  std::vector<RECORDLIVETIME> LoadLiveTime(SKSNSIMENUM::SKPERIOD);
  std::vector<RECORDLIVETIME> LoadLiveTime(int /* run-begin */,int /* run-end */);
  std::vector<RECORDLIVETIME> LoadLiveTime(std::string);

  std::tuple<int /*runnum*/, int /*num_ev*/> ConvertExpectedEvt(TRandom &, RECORDLIVETIME, double /* num_ev_per_24h */);
  std::vector<std::tuple<int /*runnum*/, int /*num_ev*/>> ConvertExpectedEvt(TRandom &, std::vector<RECORDLIVETIME>, double /* num_ev_per_24h */);

}
#endif
