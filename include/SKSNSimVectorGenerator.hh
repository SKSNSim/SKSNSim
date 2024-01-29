/*********************************
 * File: SKSNSimVectorGenerator.hh
 * Date: Wed Jul 20 23:05:48 JST 2022
 * Author: Shota Izumiyama (izumiyama@hep.phys.titech.ac.jp)
 * Desctiption:
 *********************************/

#ifndef SKSNSIMVECTORGENERATOR_H_INCLUDED
#define SKSNSIMVECTORGENERATOR_H_INCLUDED

#include <memory>
#include <TRandom3.h>
#include <algorithm>
#ifdef SKINTERNAL
#include <mcinfo.h>
#endif
#include "SKSNSimFlux.hh"
#include "SKSNSimCrosssection.hh"
#include "SKSNSimEnum.hh"
#include "SKSNSimTools.hh"

extern "C" {
  // From SKOFL
	void sn_sundir_( int *, int *, float *, float *, float *);
}

//======================
// Useful tempolary container for Linear Algebra
// These should be used in only SKSNSimVectorGenerator class internally
template<class T>
class UtilVector3 {
  public:
    T x,y,z;
    UtilVector3(const T xx, const T yy, const T zz) : x(xx), y(yy), z(zz){}
    UtilVector3(const T xyz[]) : x(xyz[0]), y(xyz[1]), z(xyz[2]){}
    UtilVector3(const T theta, const T phi) : x(std::sin(theta) * std::cos(phi)),
                                              y(std::sin(theta) * std::sin(phi)),
                                              z(std::cos(theta)){}
    T operator[](size_t i)const { if(i==0) return x; if(i==1) return y; if(i==2) return z; return 9999.;}
    T operator*(const UtilVector3 &v) const { return v.x * x + v.y * y + v.z*z; }
    UtilVector3 operator+(const UtilVector3 v) const { return UtilVector3( v.x + x,  v.y+ y, v.z+z); }
    UtilVector3 operator-(const UtilVector3 v) const { return -1.0 * v + *this; }
    T Mag2() const { return *this * *this; }
    T Mag() const { return std::sqrt(Mag2()); }
    UtilVector3 Unit() const { return (1./Mag())* *this ; }
};
template <class T> UtilVector3<T> operator*(T a, const UtilVector3<T> &v){ return UtilVector3<T>(a*v.x, a*v.y, a*v.z);}

template<class T>
class UtilMatrix3 {
  public:
    T v[3][3];
    UtilMatrix3(){}
    ~UtilMatrix3(){}
    UtilMatrix3(
        T x00, T x01, T x02,
        T x10, T x11, T x12,
        T x20, T x21, T x22
        ){
      v[0][0] = x00; v[0][1] = x01; v[0][2] = x02;
      v[1][0] = x10; v[1][1] = x11; v[1][2] = x12;
      v[2][0] = x20; v[2][1] = x21; v[2][2] = x22;
    }
    UtilMatrix3(
        UtilVector3<T> v0, UtilVector3<T> v1, UtilVector3<T> v2
        ){
      v[0][0] = v0.x; v[0][1] = v1.x; v[0][2] = v2.x;
      v[1][0] = v0.y; v[1][1] = v1.y; v[1][2] = v2.y;
      v[2][0] = v0.z; v[2][1] = v1.z; v[2][2] = v2.z;
    }
    UtilVector3<T> operator*(const UtilVector3<T> &l) const {
      T x[3];
      for(int i = 0; i < 3; i++){
        x[i] = UtilVector3<T>(v[i][0], v[i][1], v[i][2]) * l;
      }
      return UtilVector3<T>(x[0], x[1], x[2]);
    }

    UtilMatrix3 operator*(const UtilMatrix3 &m) const {
      auto r = UtilMatrix3();
      for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
          r[i][j] = 0;
          for(int k = 0; k < 3; k++)
            r[i][j] += this->v[i][k] * m.v[k][j];
        }
      }
      return r;
    }
};
//====================================

enum SKSNSimReactionType { mReactionTypeNuebarIBD = 0,
  mReactionTypeNueElastic, mReactionTypeNuebarElastic,
  mReactionTypeNuxElastic, mReactionTypeNuxbarElastic,
  mReactionTypeNueO, mReactionTypeNuebarO,
  mReactionTypeNueOsub, mReactionTypeNuebarOsub,
  mReactionTypeNuPNC, mReactionTypeNubarPNC,
  mReactionTypeNuNNC, mReactionTypeNubarNNC,
  mNReactionType };

class SKSNSimSNEventVector {
  // based on MCInfo
  private:
    int m_runnum;
    int m_subrunnum;
    
    size_t m_n_randomthrow;
    double m_weight_maxprob;
    double m_weight;

    struct VERTEX {
      double x,y,z; // cm
      int iflvvc;
      int iparvc;
      double time; // ns
    };
    std::vector<VERTEX> m_vertexs;

    struct TRACK {
      int pid; // PID code of PDG
      double energy; // total enervy  [MeV]
      double momentum_x, momentum_y, momentum_z; // MeV
      int iorgvc;
      int ivtivc;
      int ivtfvc;
      int iflgvc;
      int icrnvc;
    };
    std::vector<TRACK> m_tracks;

    struct SNINFO {
      int iEvt;
      int rType;
      double rTime;
      double rVtx[3];
      int nuType;
      double nuEne;
      double nuDir[3];
    } sninfo;

		///static void determineKinematics( SKSNSimSNEventVector &p);

    // Configuration of generator at exectuted-time
    unsigned int m_randomseed;

  public:
    SKSNSimSNEventVector() : m_n_randomthrow(0), m_weight_maxprob(0.), m_weight(0.) {sninfo.iEvt = -1;};
    ~SKSNSimSNEventVector() {};
    int AddVertex(
        double x, double y, double z,
        int iflvvc, int iparvc,
        double time){ m_vertexs.push_back(VERTEX{x,y,z,iflvvc, iparvc, time}); return m_vertexs.size();}
    int AddTrack(
        int pid, double energy, // pid, energy
        double px, double py, double pz, // momentum x,y,z
        int iorgvc, int ivtivc, int ivtfvc, int iflgvc, int icrnvc){ m_tracks.push_back(TRACK{pid,energy, px, py, pz, iorgvc, ivtivc, ivtfvc, iflgvc, icrnvc}); return m_tracks.size(); }
    int SetRunnum(int r) { m_runnum = r; return m_runnum; }
    int SetSubRunnum(int sr) { m_subrunnum = sr; return m_subrunnum; }
    int GetRunnum() const { return m_runnum; }
    int GetSubRunnum() const { return m_subrunnum; }
    int GetNVertex() const { return m_vertexs.size(); }
    double GetVertexPositionX (int i) const { if ( i >= 0 && i < GetNVertex()) return m_vertexs[i].x; return 9999.;}
    double GetVertexPositionY (int i) const { if ( i >= 0 && i < GetNVertex()) return m_vertexs[i].y; return 9999.;}
    double GetVertexPositionZ (int i) const { if ( i >= 0 && i < GetNVertex()) return m_vertexs[i].z; return 9999.;}
    int GetVertexIFLVVC (int i) const { if ( i >= 0 && i < GetNVertex()) return m_vertexs[i].iflvvc; return 9999;}
    int GetVertexIPARVC (int i) const { if ( i >= 0 && i < GetNVertex()) return m_vertexs[i].iparvc; return 9999;}
    double GetVertexTime   (int i) const { if ( i >= 0 && i < GetNVertex()) return m_vertexs[i].time; return 9999.;}
    int GetNTrack() const { return m_tracks.size(); }
    int GetTrackPID(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].pid; return 0; }
    double GetTrackEnergy(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].energy; return 0; }
    double GetTrackMomentumX(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].momentum_x; return 9999.; }
    double GetTrackMomentumY(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].momentum_y; return 9999.; }
    double GetTrackMomentumZ(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].momentum_z; return 9999.; }
    int GetTrackIORGVC(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].iorgvc; return 9999; }
    int GetTrackIVTIVC(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].ivtivc; return 9999; }
    int GetTrackIVTFVC(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].ivtfvc; return 9999; }
    int GetTrackIFLGVC(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].iflgvc; return 9999; }
    int GetTrackICRNVC(int i) const { if( i>= 0 && i < GetNTrack()) return m_tracks[i].icrnvc; return 9999; }

    void SetSNEvtInfo(const int nReact, const double tReact, const int nuType, const double nuEne, const double snDir[3], const double rVtx[3]){
      sninfo.rType = nReact;
      sninfo.rTime = tReact;
      sninfo.nuType = nuType;
      sninfo.nuEne = nuEne;
      sninfo.nuDir[0] = snDir[0];
      sninfo.nuDir[1] = snDir[1];
      sninfo.nuDir[2] = snDir[2];
      sninfo.rVtx[0] = rVtx[0];
      sninfo.rVtx[1] = rVtx[1];
      sninfo.rVtx[2] = rVtx[2];
    };
    void SetSNEvtInfoIEvt(int i) { sninfo.iEvt = i; }
    auto GetSNEvtInfoIEvt() const { return sninfo.iEvt; }
    auto GetSNEvtInfoRVtx(int i) const { return sninfo.rVtx[i]; }
    auto GetSNEvtInfoNuEne() const { return sninfo.nuEne; }
    auto GetSNEvtInfoRType() const { return sninfo.rType; }
    auto GetSNEvtInfoRTime() const { return sninfo.rTime; }
    auto GetSNEvtInfoNuType() const { return sninfo.nuType; }
    UtilVector3<double> GetSNEvtInfoNuDir() const { return UtilVector3<double>(sninfo.nuDir); }
    auto GetSNEvtInfoNuDir(int i) const { return sninfo.nuDir[i]; }

    auto GetNRandomThrow() const { return m_n_randomthrow; }
    auto SetNRandomThrow(size_t n) { m_n_randomthrow = n; return GetNRandomThrow(); }
    auto AddNRandomThrow(size_t n) { m_n_randomthrow += n; return GetNRandomThrow(); }

    auto GetWeightMaxProb() const { return m_weight_maxprob; }
    auto SetWeightMaxProb(const double w) { m_weight_maxprob = w; return GetWeightMaxProb(); }

    auto GetWeight() const { return m_weight; }
    auto SetWeight(const double w) { m_weight= w; return GetWeight(); }

    unsigned int GetRandomSeed() const { return m_randomseed; }
    unsigned int SetRandomSeed(const unsigned int r) { m_randomseed = r; return GetRandomSeed(); } // This does NOT apply seed value. Just holding runtime-configuration

    bool operator< (const SKSNSimSNEventVector &a){ return sninfo.rTime < a.sninfo.rTime; }
};

class SKSNSimVectorGenerator {
  private:
    std::vector<std::unique_ptr<SKSNSimFluxModel>> fluxmodels;
    std::vector<std::unique_ptr<SKSNSimCrosssectionModel>> xsecmodels;
    SKSNSimSNEventVector GenerateSNEvent(){
      return SKSNSimSNEventVector();
    }

    //=====================
    // Configuration
    double m_generator_energy_min;
    double m_generator_energy_max;
    bool m_runtime_norm;
    double m_runtime_factor;
    int m_runtime_runbegin;
    int m_runtime_runend;
    int m_runtime_period;
    int m_runnum;
    int m_subrunnum;
    bool m_flat_pos_energy;
    SKSNSIMENUM::TANKVOLUME m_generator_volume;
    std::shared_ptr<TRandom> randomgenerator;
    unsigned int m_randomseed;
    //===================== end configuration


    // For hit-and-miss method
    double m_max_hit_probability; // maximum of (flux) x (xsec) // should be updated with new flux or xsec models
    static double FindMaxProb ( SKSNSimFluxModel &, SKSNSimCrosssectionModel &, int /* elapseday */ = -1);
    double SetMaximumHitProbability();
    
  public:
    SKSNSimVectorGenerator():
      m_runnum((int)SKSNSIMENUM::SKPERIODRUN::SKMC),
      m_subrunnum(0),
      m_flat_pos_energy ( false ),
      m_generator_volume(SKSNSIMENUM::TANKVOLUME::kIDFULL)
    {}
    ~SKSNSimVectorGenerator(){}
    void AddFluxModel(SKSNSimFluxModel *fm){ fluxmodels.push_back(std::move(std::unique_ptr<SKSNSimFluxModel>(fm))); SetMaximumHitProbability(); } // after this, the pointer will be managed by SKSNSimVectorGenerator class
    void AddXSecModel(SKSNSimCrosssectionModel *xm){ xsecmodels.push_back(std::move(std::unique_ptr<SKSNSimCrosssectionModel>(xm))); SetMaximumHitProbability(); } // after this, the pointer will be managed by SKSNSimVectorGenerator class
    SKSNSimSNEventVector GenerateEventIBD();
    SKSNSimSNEventVector GenerateEventIBDFlat();
    SKSNSimSNEventVector GenerateEvent() { return m_flat_pos_energy? GenerateEventIBDFlat(): GenerateEventIBD(); }; // Tentatively, supporting only IBD channel
    std::vector<SKSNSimSNEventVector> GenerateEvents(int n) {
      std::vector<SKSNSimSNEventVector> buf(n);
      for(auto it = buf.begin(); it != buf.end(); it++) *it = GenerateEvent();
      return buf;
    } ;
    // std::vector<SKSNSimSNEventVector> GenerateEventsAlongLivetime(int, int, double);

    // Configuration
    double SetEnergyMin(const double e){ m_generator_energy_min = e; return m_generator_energy_min;}
    double SetEnergyMax(const double e){ m_generator_energy_max = e; return m_generator_energy_max;}
    double GetEnergyMin() const {return m_generator_energy_min;}
    double GetEnergyMax() const {return m_generator_energy_max;}
    SKSNSIMENUM::TANKVOLUME GetGeneratorVolume() const {return m_generator_volume;}
    SKSNSIMENUM::TANKVOLUME SetGeneratorVolume(const SKSNSIMENUM::TANKVOLUME v){ m_generator_volume = v; return GetGeneratorVolume();}
    bool SetNormRuntime(const bool f) { m_runtime_norm = f; return m_runtime_norm; }
    bool GetNormRuntime() const { return m_runtime_norm; }
    double  SetRuntimeFactor(const double r) { m_runtime_factor = r; return m_runtime_factor; }
    double  GetRuntimeFactor() const { return m_runtime_factor; }
    int SetRuntimeBegin(const int r) { m_runtime_runbegin = r; return m_runtime_runbegin; }
    int GetRuntimeBegin() const { return m_runtime_runbegin; }
    int SetRuntimeEnd(const int r) { m_runtime_runend = r; return m_runtime_runend; }
    int GetRuntimeEnd() const { return m_runtime_runend; }
    int SetRuntimePeriod(const int r) { m_runtime_period = r; return m_runtime_period; }
    int GetRuntimePeriod() const { return m_runtime_period; }
    int SetRUNNUM(const int r) { m_runnum = r; return m_runnum; }
    int GetRUNNUM() const { return m_runnum; }
    int SetSubRUNNUM(const int r) { m_subrunnum = r; return m_subrunnum; }
    int GetSubRUNNUM() const { return m_subrunnum; }
    unsigned int GetRandomSeed() const {return m_randomseed; }
    unsigned int SetRandomSeed(unsigned int s) { m_randomseed = s; return GetRandomSeed(); } // This does NOT apply the seed. Just holding the runtime-information.
    void   SetRandomGenerator(std::shared_ptr<TRandom> rng) { randomgenerator = rng; }
    bool SetFlatPositronFlux(const bool f) { m_flat_pos_energy = f; return m_flat_pos_energy; }
    bool GetFlatPositronFlux() const { return m_flat_pos_energy; }
};

class SKSNSimVectorSNGenerator {
  private:
    std::vector<std::unique_ptr<SKSNSimFluxModel>> fluxmodels;
    std::map<XSECTYPE, std::shared_ptr<SKSNSimCrosssectionModel>> xsecmodels;
    SKSNSimSNEventVector GenerateSNEvent(){
      return SKSNSimSNEventVector();
    }
    // Event header
    int m_runnum;
    int m_subrunnum;

    // Vector generator
    double m_generator_energy_min;
    double m_generator_energy_max;
    size_t m_nu_energy_nbins;
    double m_generator_time_min;
    double m_generator_time_max;
    size_t m_time_nbins;
    bool   m_fill_event;
    SKSNSIMENUM::TANKVOLUME m_generator_volume;

    // Phsics assumption
    SKSNSIMENUM::NEUTRINOOSCILLATION m_nuosc_type;

    // SN direction
    int m_sn_date[3];
    int m_sn_time[3];
    double m_sn_dir[3];
    double m_distance_kpc;

    // double m_max_hit_probability; // maximum of (flux) x (xsec) // should be updated with new flux or xsec models

    std::shared_ptr<TRandom> randomgenerator;
    unsigned int m_randomseed;

    static double FindMaxProb ( const double, const SKSNSimCrosssectionModel &);

    //double SetMaximumHitProbability();
    std::vector<SKSNSimSNEventVector> MakeEvent(const double nuEneBinSize, const double tBinSize, const double time, const double nu_energy, const int nReact, const int nuType, const double rate);
    void FillEvent(std::vector<SKSNSimSNEventVector> &evt_buffer);
    static void determineKinematics( std::map<XSECTYPE, std::shared_ptr<SKSNSimCrosssectionModel>> xsecmodels, TRandom &rng, SKSNSimSNEventVector &ev, const double snDir[]);
    static void determineKinematicsIBD( const SKSNSimXSecIBDSV & xsec, TRandom &rng, SKSNSimSNEventVector &ev, const UtilVector3<double> nuDir);
    static void determineAngleNuebarP( TRandom &rng, const SKSNSimXSecIBDSV & xsec, const double nuEne, double & eEne, double & eTheta, double & ePhi );
    static void determineAngleElastic( TRandom &rng, const SKSNSimXSecNuElastic & xsec, const int nReact, const double nuEne, double & eEne, double & eTheta, double & ePhi, int &iSkip );
    static void determineAngleNueO(TRandom &rng, SKSNSimXSecNuOxygen &xsec, const int Reaction, const int State, const int Ex_state, const int channel, const double nuEne, double & eEne, double & eTheta, double & ePhi );

    
  public:
    SKSNSimVectorSNGenerator();
    ~SKSNSimVectorSNGenerator(){ SKSNSimTools::DumpDebugMessage(" dtor of SKSNSimVectorSNGenerator");}
    void AddFluxModel(SKSNSimFluxModel *fm){ fluxmodels.push_back(std::move(std::unique_ptr<SKSNSimFluxModel>(fm))); /* SetMaximumHitProbability(); */ } // after this, the pointer will be managed by SKSNSimVectorGenerator class
    void AddFluxModel(std::unique_ptr<SKSNSimFluxModel> fm){ fluxmodels.push_back(std::move(fm)); /* SetMaximumHitProbability(); */ } // after this, the pointer will be managed by SKSNSimVectorGenerator class
    std::vector<SKSNSimSNEventVector> GenerateEvents();

    //========================================
    // Configuration
    double SetEnergyMin(const double e){ m_generator_energy_min = e; return m_generator_energy_min;}
    double SetEnergyMax(const double e){ m_generator_energy_max = e; return m_generator_energy_max;}
    size_t SetEnergyNBins(const size_t n) { m_nu_energy_nbins = n; return GetEnergyNBins(); }
    double GetEnergyMin() const {return m_generator_energy_min;}
    double GetEnergyMax() const {return m_generator_energy_max;}
    size_t GetEnergyNBins() const { return m_nu_energy_nbins; } 
    double GetEnergyBinWidth() const { return (GetEnergyMax() - GetEnergyMin())/(double)GetEnergyNBins(); } 
    double SetTimeMin(const double e){ m_generator_time_min = e; return m_generator_time_min;}
    double SetTimeMax(const double e){ m_generator_time_max = e; return m_generator_time_max;}
    size_t SetTimeNBins(const size_t n) { m_nu_energy_nbins = n; return GetTimeNBins(); }
    double GetTimeMin() const {return m_generator_time_min;}
    double GetTimeMax() const {return m_generator_time_max;}
    size_t GetTimeNBins() const {return m_time_nbins;}
    double GetTimeBinWidth() const { return (GetTimeMax() - GetTimeMin())/(double)GetTimeNBins(); }
    bool   GetFlagFillEvent() const { return m_fill_event; }
    bool   SetFlagFillEvent(const bool f){ m_fill_event = f; return GetFlagFillEvent(); }
    unsigned int GetRandomSeed() const {return m_randomseed; }
    unsigned int SetRandomSeed(unsigned int s) { m_randomseed = s; return GetRandomSeed(); } // This does NOT apply the seed. Just holding the runtime-information.
    void   SetRandomGenerator(std::shared_ptr<TRandom> rng) { randomgenerator = rng; }
    SKSNSIMENUM::TANKVOLUME GetGeneratorVolume() const { return m_generator_volume; }
    SKSNSIMENUM::TANKVOLUME SetGeneratorVolume(SKSNSIMENUM::TANKVOLUME v) { m_generator_volume = v; return GetGeneratorVolume(); }
    SKSNSIMENUM::TANKVOLUME SetGeneratorVolume(int v) { m_generator_volume = (SKSNSIMENUM::TANKVOLUME)v; return GetGeneratorVolume(); }
    double GetSNDistanceKpc() const { return m_distance_kpc;}
    double GetSNDistanceRatioTo10kpc() const { return  pow(10.0 / GetSNDistanceKpc(),2.); }
    double SetSNDistanceKpc(const double d) { m_distance_kpc = d; return GetSNDistanceKpc();}
    SKSNSIMENUM::NEUTRINOOSCILLATION GetGeneratorNuOscType () const { return m_nuosc_type; }
    SKSNSIMENUM::NEUTRINOOSCILLATION SetGeneratorNuOscType(SKSNSIMENUM::NEUTRINOOSCILLATION t) { m_nuosc_type = t; return GetGeneratorNuOscType(); }
    SKSNSIMENUM::NEUTRINOOSCILLATION SetGeneratorNuOscType(int t) { m_nuosc_type = (SKSNSIMENUM::NEUTRINOOSCILLATION)t; return GetGeneratorNuOscType(); }
    int GetRUNNUM() const { return m_runnum; }
    int SetRUNNUM(const int r){ m_runnum = r; return GetRUNNUM(); }
    int GetSubRUNNUM() const { return m_subrunnum; }
    int SetSubRUNNUM(const int r){ m_subrunnum = r; return GetSubRUNNUM(); }
};

#endif
