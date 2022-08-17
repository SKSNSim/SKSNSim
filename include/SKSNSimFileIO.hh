/**********************************
 * File: SKSNSimFileIO.hh
 **********************************/

#ifndef __SKSNSIMFILEIO_H_INCLUDED__
#define __SKSNSIMFILEIO_H_INCLUDED__

#include <TFile.h>
#include <TTree.h>
#include <mcinfo.h>
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
    virtual void Write(const std::vector<SKSNSimSNEventVector> &vecs) {for(auto it = vecs.begin(); it != vecs.end(); it++) Write(*it);};
};

std::vector<SKSNSimFileSet> GenerateOutputFileList(SKSNSimUserConfiguration &conf);

class SKSNSimFileOutTFile : SKSNSimFileOutput {
  private:
    TFile *m_fileptr;
    MCInfo *m_MC;
    TTree *m_OutTree;


  public:
    SKSNSimFileOutTFile (){ m_fileptr = NULL; }
    SKSNSimFileOutTFile (const std::string fname) { Open(fname); }
    ~SKSNSimFileOutTFile (){ if (m_fileptr != NULL) delete m_fileptr; }

    void Open(const std::string fname);
    void Close();

    void Write(const SKSNSimSNEventVector &);
};
#endif
