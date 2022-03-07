#include <limits.h>
#include <string>
#include <TApplication.h>
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <stdarg.h>
#include <stdlib.h>
#include <time.h>
#include <TROOT.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TMinuit.h>
#include <TNtuple.h>
#include <TNetFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF2.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TStyle.h>
#include "TTree.h"
#include "TChain.h" 
#include <TTreeCache.h>

void event_dist()
{

  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0);

  //TCanvas *c0 = new TCanvas("c0","c0",10,10,750,900);
  TCanvas *c0 = new TCanvas("c0","c0",10,10,750,750);
  c0->SetGrid();
  c0->SetLogx();
  c0->SetLogy();
  c0->SetLeftMargin(0.2);
  c0->SetBottomMargin(0.13);
  c0->SetRightMargin(0.1);
  c0->SetTopMargin(0.2);
  c0->SetTicks(1,1);
  c0->SetFrameLineWidth(3);

  TH1F *frame1 = gPad->DrawFrame(0.1,0.03,5000,300000000);
  frame1->GetXaxis()->SetTitleSize(0.05);
  frame1->GetXaxis()->SetLabelSize(0.05);
  frame1->GetXaxis()->SetTitleOffset(1.3);
  frame1->GetXaxis()->SetTitle("distance (kpc)");

  frame1->GetYaxis()->SetTitleSize(0.05);
  frame1->GetYaxis()->SetLabelSize(0.05);
  frame1->GetYaxis()->SetTitleOffset(1.3);
  frame1->GetYaxis()->SetTitle("events/32.5 kton");

  double totNuebarp, totNueElastic, totNuebarElastic, totNuxElastic, totNuxbarElastic;
  double totNum;
  TF1 *f1;

  // Nakazato models

  string str_arr[22] = {"../VectGen/data/intp1301.data.dat", "../VectGen/data/intp1302.data.dat","../VectGen/data/intp1303.data.dat","../VectGen/data/intp1311.data.dat","../VectGen/data/intp1312.data.dat","../VectGen/data/intp1313.data.dat","../VectGen/data/intp2001.data.dat","../VectGen/data/intp2002.data.dat","../VectGen/data/intp2003.data.dat","../VectGen/data/intp2011.data.dat","../VectGen/data/intp2012.data.dat","../VectGen/data/intp2013.data.dat","../VectGen/data/intp3001.data.dat","../VectGen/data/intp3002.data.dat","../VectGen/data/intp3003.data.dat","../VectGen/data/intp5001.data.dat","../VectGen/data/intp5002.data.dat","../VectGen/data/intp5003.data.dat","../VectGen/data/intp5011.data.dat","../VectGen/data/intp5012.data.dat","../VectGen/data/intp5013.data.dat","../VectGen/data/spectob3010.data.dat"};

  int maxModel, minModel;
  double maxNum=0., minNum=999999., disBete = 0.15/10., disLMC = 50./10., disM31 = 740./10.;
  for(int i=0; i<22;i++){
    ifstream data(str_arr[i].c_str());
    data >> totNuebarp >> totNueElastic >> totNuebarElastic >> totNuxElastic >> totNuxbarElastic;
    totNum = totNuebarp + totNueElastic + totNuebarElastic + totNuxElastic + totNuxbarElastic;
    std::cout << i << " " << totNum << std::endl;
    if(i<21) {
      if(totNum > maxNum) {maxNum=totNum; maXModel=i;}
      if(totNum < minNum) {minNum=totNum; minModel=i;}
    }
    else {
      std::cout << totNum / (disBete * disBete) << " " <<totNum / (disLMC * disLMC) << " " <<totNum / (disM31 * disM31) << std::endl;
    }
  
    f1 = new TF1("f1","[0]/x/x",0.01, 5000);
    f1->SetParameter(0, totNum*100.);
    f1->SetLineWidth(0.05);
    f1->Draw("same");
  }
  std::cout << "Minumum" << std::endl;
  std::cout << minModel << " " << minNum << std::endl;
  std::cout << minNum / (disBete * disBete) << " " <<minNum / (disLMC * disLMC) << " " <<minNum / (disM31 * disM31) << std::endl;
  std::cout << "Maximum" << std::endl;
  std::cout << maxModel << " " << maxNum << std::endl;
  std::cout << maxNum / (disBete * disBete) << " " <<maxNum / (disLMC * disLMC) << " " <<maxNum / (disM31 * disM31) << std::endl;


  // Input manually for Livermore (Totani)

  double NumLiv = 9398.64;
  std::cout << "Livermore" << std::endl;
  std::cout << NumLiv / (disBete * disBete) << " " << NumLiv << " " <<NumLiv / (disLMC * disLMC) << " " <<NumLiv / (disM31 * disM31) << std::endl;

  f1 = new TF1("f1","[0]/x/x",0.01, 5000);
  f1->SetParameter(0, NumLiv*100.);
  f1->SetLineWidth(0.05);
  f1->SetLineColor(2);
  f1->Draw("same");

  // Write the candidate sources

  //double xx[4]={0.15, 8., 50., 740.}, yy[4]={270000000.,270000000.,270000000.,270000000.};
  double xx[4]={0.15, 10., 50., 740.}, yy[4]={270000000.,270000000.,270000000.,270000000.};
  TGraph *gr1 = new TGraph(4,xx,yy);
  gr1->SetMarkerStyle(29);
  gr1->SetMarkerSize(2);
  gr1->Draw("same,p");

  TText *t;
  t = new TText(0.15, 380000000, "Betelgeuse");
  t->SetTextAlign(2);
  t->SetTextAngle(90);
  t->SetTextSize(0.03);
  t->Draw("same");

  t = new TText(0.2, 380000000, "Antares");
  t->SetTextAlign(2);
  t->SetTextAngle(90);
  t->SetTextSize(0.03);
  t->Draw("same");

  t = new TText(8., 380000000, "Galactic");
  t->SetTextAlign(2);
  t->SetTextAngle(90);
  t->SetTextSize(0.03);
  t->Draw("same");

  t = new TText(11., 380000000, "center");
  t->SetTextAlign(2);
  t->SetTextAngle(90);
  t->SetTextSize(0.03);
  t->Draw("same");

  t = new TText(50., 380000000, "LMC");
  t->SetTextAlign(2);
  t->SetTextAngle(90);
  t->SetTextSize(0.03);
  t->Draw("same");

  t = new TText(740., 380000000, "M31");
  t->SetTextAlign(2);
  t->SetTextAngle(90);
  t->SetTextSize(0.03);
  t->Draw("same");

  c0->Print("event_dist.pdf");

}
