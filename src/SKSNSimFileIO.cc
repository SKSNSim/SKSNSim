/**********************************
 * File: SKSNSimFileIO.cc
 **********************************/

#include "SKSNSimFileIO.hh"
#include <TFile.h>
#include <TFileCacheWrite.h>
#include <TTree.h>
#include <mcinfo.h>


void SKSNSimFileOutTFile::Open(const std::string fname) {
  if(m_fileptr != NULL){
    Close();
    delete m_fileptr;
  }

  m_fileptr = new TFile(fname.c_str(), "RECREATE");
	//------------------------------------------------------------------------
	// set write cache to 40MB 
	TFileCacheWrite *cachew = new TFileCacheWrite(m_fileptr, 40*1024*1024);
	//------------------------------------------------------------------------

	/*-----define class-----*/
	m_MC = new MCInfo;
	m_MC->Clear();
	m_MC->SetName("MC");

	/*-----define branch-----*/
	TList *TopBranch = new TList;
	TopBranch->Add(m_MC);

	// define tree
	m_OutTree = new TTree("data", "SK 5 tree");
	// new MF
	m_OutTree->SetCacheSize(40*1024*1024);
	int bufsize = 8*1024*1024;      // may be this is the best 15-OCT-2007 Y.T.
	m_OutTree->Branch(TopBranch,bufsize);
}

void SKSNSimFileOutTFile::Close(){
	m_fileptr->cd();
	m_OutTree->Write();
  m_fileptr->Save();
	m_OutTree->Reset();
	m_fileptr->Close();
}

void SKSNSimFileOutTFile::Write(const SKSNSimSNEventVector &ev){

  m_MC->mcinfo[0] = ev.GetRunnum();
  m_MC->mcinfo[1] = ev.GetSubRunnum();

  // MCVERTEX (see $SKOFL_ROOT/inc/vcvrtx.h )
  m_MC->nvtxvc = ev.GetNVertex();
  for(int i = 0; i < ev.GetNVertex(); i++){
    m_MC->pvtxvc[i][0] = ev.GetVertexPositionX(i);
    m_MC->pvtxvc[i][1] = ev.GetVertexPositionY(i);
    m_MC->pvtxvc[i][2] = ev.GetVertexPositionZ(i);
    m_MC->timvvc[i] = ev.GetVertexTime(i);
    m_MC->iflvvc[i] = ev.GetVertexIFLVVC(i);
    m_MC->iparvc[i] = ev.GetVertexIPARVC(i);
  }


  m_MC->nvc = ev.GetNTrack();
  for(int i = 0; i < ev.GetNTrack(); i++){
    m_MC->ipvc[i] = ev.GetTrackPID(i);
    m_MC->energy[i] = ev.GetTrackEnergy(i);
    m_MC->pvc[i][0] = ev.GetTrackMomentumX(i);
    m_MC->pvc[i][1] = ev.GetTrackMomentumY(i);
    m_MC->pvc[i][2] = ev.GetTrackMomentumZ(i);
    m_MC->iorgvc[i] = ev.GetTrackIORGVC(i);
    m_MC->ivtivc[i] = ev.GetTrackIVTIVC(i);
    m_MC->ivtfvc[i] = ev.GetTrackIVTFVC(i);
    m_MC->iflgvc[i] = ev.GetTrackIFLGVC(i);
    m_MC->icrnvc[i] = ev.GetTrackICRNVC(i);
  }

  m_OutTree->Fill();
}
