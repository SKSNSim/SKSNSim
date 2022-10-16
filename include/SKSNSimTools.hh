#ifndef __SKSNSIMTOOLS_H_INCLUDED__
#define __SKSNSIMTOOLS_H_INCLUDED__

#include <TString.h>
#include <iostream>

namespace SKSNSimTools{
  auto DumpDebugMessage = [] (TString str) {
#ifdef DEBUG
    std::cout << "[IZU:DEBUG] " << str << std::endl;
#endif
  };
}
#endif
