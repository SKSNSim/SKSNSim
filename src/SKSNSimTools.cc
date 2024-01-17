
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "SKSNSimTools.hh"

extern "C" {
  void elapseday_date_(int *, int *, int *, int *);
  void elapseday_run_(int *, int *);
}

namespace SKSNSimTools {
  SKSNSIMENUM::SKPERIOD FindSKPeriod(int rn /* run_number */) {
    auto checkRange = [] (int t, int b, int e) {
      return ( t >= b && t < e );
    };
    if( checkRange(rn, (int)SKSNSIMENUM::SKPERIODRUN::SKIBEGIN, (int)SKSNSIMENUM::SKPERIODRUN::SKIEND) )
      return SKSNSIMENUM::SKPERIOD::SKI;
    else if( checkRange(rn, (int)SKSNSIMENUM::SKPERIODRUN::SKIIBEGIN, (int)SKSNSIMENUM::SKPERIODRUN::SKIIEND) )
      return SKSNSIMENUM::SKPERIOD::SKII;
    else if( checkRange(rn, (int)SKSNSIMENUM::SKPERIODRUN::SKIIIBEGIN, (int)SKSNSIMENUM::SKPERIODRUN::SKIIIEND) )
      return SKSNSIMENUM::SKPERIOD::SKIII;
    else if( checkRange(rn, (int)SKSNSIMENUM::SKPERIODRUN::SKIVBEGIN, (int)SKSNSIMENUM::SKPERIODRUN::SKIVEND) )
      return SKSNSIMENUM::SKPERIOD::SKIV;
    else if( checkRange(rn, (int)SKSNSIMENUM::SKPERIODRUN::SKVBEGIN, (int)SKSNSIMENUM::SKPERIODRUN::SKVEND) )
      return SKSNSIMENUM::SKPERIOD::SKV;
    else if( checkRange(rn, (int)SKSNSIMENUM::SKPERIODRUN::SKVIBEGIN, (int)SKSNSIMENUM::SKPERIODRUN::SKVIEND) )
      return SKSNSIMENUM::SKPERIOD::SKVI;
    else if( checkRange(rn, (int)SKSNSIMENUM::SKPERIODRUN::SKVIIBEGIN, (int)SKSNSIMENUM::SKPERIODRUN::SKVIIEND) )
      return SKSNSIMENUM::SKPERIOD::SKVII;

    /* default */
    return SKSNSIMENUM::SKPERIOD::NSKPERIOD;
  }

  int elapseday(int yy, int mm, int dd){
    int eladay = -1;
    elapseday_date_(&yy, &mm, &dd, &eladay);
    return eladay;
  }
  int elapseday(int run){
    int eladay = -1;
    elapseday_run_(&run, &eladay);
    return eladay;
  }
}

namespace SKSNSimLiveTime {
  const static std::map<SKSNSIMENUM::SKPERIOD, std::string> FNAMEMAP = {
    { SKSNSIMENUM::SKPERIOD::SKV, "/home/sklowe/realtime_sk5_rep/solar_oct19/livetime/livetime5.r080539.r082086.txt" },
    { SKSNSIMENUM::SKPERIOD::SKVI, "/home/sklowe/realtime_sk6_rep/solar_oct22/livetime/livetime5.r085220.r087220.txt" }
    { SKSNSIMENUM::SKPERIOD::SKVII, "/home/sklowe/realtime_sk7_rep/solar_nov23/livetime/livetime5.r080000.r091985.txt" }
  };

  std::vector<std::tuple<int,double>> LoadLiveTime(std::string fname){
    std::vector<std::tuple<int, double>> buf;

    std::cout << "LoadLiveTime: open livetime file: " << fname << std::endl;
    std::ifstream flive ( fname );
    if( flive.is_open() ){
      for( std::array<char, 500> a; flive.getline(&a[0], 500); ){
        std::stringstream ss ( std::string(a.data(), 500) );
        double run, liveall, liveday, livenight;
        ss >> run >> liveall >> liveday >> livenight;
        if(run >= 0 && run < 999999) buf.push_back( std::tuple<int,double>(run, liveall));
      }
    }
    flive.close();

    return buf;
  }

  std::vector<std::tuple<int,double>> LoadLiveTime(SKSNSIMENUM::SKPERIOD p){
    if( FNAMEMAP.find(p) == FNAMEMAP.end() ){
      std::cout << "LoadLiveTime: not supported this period: " << (int)p << std::endl;
      return std::vector<std::tuple<int, double>>();
    }

    const std::string fname = FNAMEMAP.at(p);

    std::vector<std::tuple<int, double>> buf = LoadLiveTime(fname);

    return buf;
  }

  std::vector<std::tuple<int,double>> LoadLiveTime(int rb, int re){
    std::vector<std::tuple<int, double>> buf;
    auto periodBegin = SKSNSimTools::FindSKPeriod(rb);
    auto periodEnd   = SKSNSimTools::FindSKPeriod(re);

    if( periodBegin > periodEnd ) return buf;

    for( int i = (int)periodBegin; i <= (int)periodEnd; i++){
      auto nb =  LoadLiveTime( (SKSNSIMENUM::SKPERIOD)i );
      buf.insert(buf.begin(),nb.begin(), nb.end());
    }

    auto runRangeCheck = [rb, re](RECORDLIVETIME r){ return (std::get<0>(r)<rb || std::get<0>(r)>=re);};

    auto removedEnd = std::remove_if(buf.begin(), buf.end(), runRangeCheck);
    buf.erase(removedEnd, buf.end());
    std::sort( buf.begin(), buf.end());
    return buf;
  }

  std::tuple<int /*runnum*/, int /*num_ev*/> ConvertExpectedEvt(TRandom &rng, RECORDLIVETIME rec, double nev_24h){
    double liveday = std::get<1>(rec);
    double evt_double = liveday * nev_24h;
    int nev = (int)evt_double; // getting integer parts
    {
      double frac = rng.Uniform();
      if( frac <= evt_double - (double)nev ) nev++;
    }

    return std::make_tuple( std::get<0>(rec), nev);
  }

  std::vector<std::tuple<int , int>> ConvertExpectedEvt(TRandom &rng, std::vector<RECORDLIVETIME> recs, double nev_24h){
    std::vector<std::tuple<int,int>> buf(recs.size());

    for(auto it = recs.begin(); it != recs.end(); it++){
      auto a = ConvertExpectedEvt(rng, *it, nev_24h);
      buf[ std::distance(recs.begin(),it) ] = a;
    }

    return buf;
  }
}
