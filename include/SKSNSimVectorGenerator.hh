/*********************************
 * File: SKSNSimVectorGenerator.hh
 * Date: Wed Jul 20 23:05:48 JST 2022
 * Author: Shota Izumiyama (izumiyama@hep.phys.titech.ac.jp)
 * Desctiption:
 *********************************/

#ifndef __SKSNSIMVECTORGENERATOR_H_INCLUDED__
#define __SKSNSIMVECTORGENERATOR_H_INCLUDED__

#include <memory>
#include <mcinfo.h>
#include <TRandom3.h>
#include "SKSNSimFlux.hh"
#include "SKSNSimCrosssection.hh"

class SKSNSimSNEventVector {
  // based on MCInfo
  private:
    int m_runnum;
    int m_subrunnum;

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


  public:
    SKSNSimSNEventVector() {};
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
};

class SKSNSimVectorGenerator {
  private:
    std::vector<std::unique_ptr<SKSNSimFluxModel>> fluxmodels;
    std::vector<std::unique_ptr<SKSNSimCrosssectionModel>> xsecmodels;
    SKSNSimSNEventVector GenerateSNEvent(){
      return SKSNSimSNEventVector();
    }

    double m_generator_energy_min;
    double m_generator_energy_max;

    double m_max_hit_probability; // maximum of (flux) x (xsec) // should be updated with new flux or xsec models

    TRandom3 randomgenerator;

    static double FindMaxProb ( SKSNSimFluxModel &, SKSNSimCrosssectionModel &);

    double SetMaximumHitProbability();
    
  public:
    SKSNSimVectorGenerator(){}
    ~SKSNSimVectorGenerator(){}
    void AddFluxModel(SKSNSimFluxModel *fm){ fluxmodels.push_back(std::move(std::unique_ptr<SKSNSimFluxModel>(fm))); SetMaximumHitProbability(); } // after this, the pointer will be managed by SKSNSimVectorGenerator class
    void AddXSecModel(SKSNSimCrosssectionModel *xm){ xsecmodels.push_back(std::move(std::unique_ptr<SKSNSimCrosssectionModel>(xm))); SetMaximumHitProbability(); } // after this, the pointer will be managed by SKSNSimVectorGenerator class
    SKSNSimSNEventVector GenerateEvent();
    SKSNSimSNEventVector GenerateEventIBD();
    std::vector<SKSNSimSNEventVector> GenerateEvents(int);
    std::vector<SKSNSimSNEventVector> GenerateEventsAlongLivetime(int, int);
    double SetEnergyMin(const double e){ m_generator_energy_min = e; return m_generator_energy_min;}
    double SetEnergyMax(const double e){ m_generator_energy_max = e; return m_generator_energy_max;}
    double GetEnergyMin() const {return m_generator_energy_min;}
    double GetEnergyMax() const {return m_generator_energy_max;}
    void SetRandomSeed(const unsigned s) { randomgenerator.SetSeed(s); }
    unsigned GetRandomSeed(const int s) { return randomgenerator.GetSeed(); }
};

#endif
