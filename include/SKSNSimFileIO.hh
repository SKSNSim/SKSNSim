/**********************************
 * File: SKSNSimFileIO.hh
 **********************************/

#ifndef SKSNSIMFILEIO_H_INCLUDED
#define SKSNSIMFILEIO_H_INCLUDED

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#ifdef SKINTERNAL
#include <mcinfo.h>
#include <snevtinfo.h>
#endif
#include "SKSNSimVectorGenerator.hh"
#include "SKSNSimUserConfiguration.hh"

class SKSNSimFileSet {
    std::string fname;
    int run;
    int subrun;
    int nev;
  public:
    SKSNSimFileSet(std::string fn = "", int r =999999, int sb = -1, int n = 0):fname(fn), run(r), subrun(sb), nev(n) {}
    ~SKSNSimFileSet() {}
    std::string GetFileName() const { return fname;}
    int GetRun() const { return run;}
    int GetSubrun() const { return subrun;}
    int GetNumEvents() const { return nev;}
};

class SKSNSimFileOutput{
  public:
    virtual ~SKSNSimFileOutput(){}

    virtual void Open(const std::string) = 0;
    virtual void Close() = 0;

    virtual void Write(const SKSNSimSNEventVector &) = 0;
    virtual void Write(const std::vector<SKSNSimSNEventVector> &vecs) = 0;
};

std::vector<SKSNSimFileSet> GenerateOutputFileList(SKSNSimUserConfiguration &conf);

class SKSNSimFileOutTFile : public SKSNSimFileOutput {
  private:
    TFile *m_fileptr;
#ifdef SKINTERNAL
    MCInfo *m_MC;
    SNEvtInfo *m_SN;
#else
    void *m_MC;
    void *m_SN;
#endif
    TTree *m_OutTree;

    Double_t weight;
    TTree *m_OutWeightTree;

  public:
    SKSNSimFileOutTFile () : m_fileptr(NULL), m_MC(NULL), m_SN(NULL) { }
    SKSNSimFileOutTFile (const std::string fname) : m_fileptr(NULL), m_MC(NULL), m_SN(NULL) { Open(fname); }
    ~SKSNSimFileOutTFile (){ if (m_fileptr != NULL) delete m_fileptr; }

    void Open(const std::string fname, const bool including_snevtinfo);
    void Open(const std::string fname) { Open(fname, false);};
    void Close();

    void Write(const SKSNSimSNEventVector &);
    void Write(const std::vector<SKSNSimSNEventVector> &vecs)
    {for(auto it = vecs.begin(); it != vecs.end(); it++) Write(*it);};
};

class SKSNSimFileOutNuance : public SKSNSimFileOutput {
  private:

    Double_t weight;
    std::unique_ptr<std::ofstream> ofs;

  public:
    SKSNSimFileOutNuance () {};
    SKSNSimFileOutNuance (const std::string fname) { Open(fname); }
    ~SKSNSimFileOutNuance (){ if( ofs->is_open()) Close(); }

    void Open(const std::string fname);
    void Close();

    void Write(const SKSNSimSNEventVector &);
    void Write(const std::vector<SKSNSimSNEventVector> &vecs)
    {for(auto it = vecs.begin(); it != vecs.end(); it++) Write(*it);};
};
#endif
