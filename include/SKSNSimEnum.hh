/************************************
 * File: SKSNSimEnum.hh
 * Author: Shota Izumiyama (izumiyama@hep.phys.titech.ac.jp)
 * Date: Fri Oct  7 10:08:11 JST 2022
 * Desctiption:
 *   Collection of enumeration for SKSNSim
 ***********************************/
#ifndef SKSNSIMENUM_H_INCLUDED
#define SKSNSIMENUM_H_INCLUDED

namespace SKSNSIMENUM {
  enum struct NEUTRINOOSCILLATION { kNONE = 0, kNORMAL, kINVERTED, kNNEUTRINOOSCILLATION};
  enum struct TANKVOLUME { kIDFV = 0, kIDFULL, kTANKFULL, kNTANKVOLUME};
  enum struct SKPERIODRUN { // PERIOD >= __BEGIN && PERIOD < __END (END means it is excluded)
    SKIBEGIN, SKIEND,
    SKIIBEGIN, SKIIEND,
    SKIIIBEGIN, SKIIIEND,
    SKIVBEGIN = 60000, SKIVEND = 80000, 
    SKVBEGIN = 80000, SKVEND = 85000,
    SKVIBEGIN = 85000, SKVIEND = 87999,
    SKVIIBEGIN = 88000, SKVIIEND = 92999,
    SKVIIIBEGIN = 93000, SKVIIIEND = 999998,
    SKMC = 999999 };
  enum struct SKPERIOD {
    SKI = 0, SKII, SKIII, SKIV, SKV, SKVI, SKVII, SKVIII, NSKPERIOD
  };
  enum struct NUINTERACTIONTYPE {
    kNUEBIBD =0,
    kNUEELASTIC, kNUEBELASTIC, kNUXELASTIC, kNUXBELASTIC,
    kNUEOXY, kNUEBOXY, kNUEOXYSUB, kNUEBOXYSUB,
    kNUNCP, kNUNCN, kNUBNCP, kNUBNCN,
    kNNUINTERACTIONTYPE
  };
};

#endif
