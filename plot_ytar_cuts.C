#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCutG.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
using namespace std;


void plot_ytar_cuts(Int_t nrun=1814,Int_t FileID=-2) {
  
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.14);
  const Int_t Number=3;
  Double_t Red[Number] = { 1.0,0.0,0.0};
  Double_t Blue[Number] = { 1.0,0.0,1.0};
  Double_t Green[Number] = { 0.0,1.0,0.0};
  Double_t Len[Number] = { 0.0,.5,1.0};
  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Len,Red,Green,Blue,nb);
  //  Get info for that optics run
  TString OpticsFile = "list_of_optics_run.dat";
  ifstream file_optics(OpticsFile.Data());
  TString opticsline;
  TString OpticsID="";
  Int_t RunNum=0.;
  Double_t CentAngle=0.;
  Int_t SieveFlag=1;
  Int_t NumFoil=0;
  TString temp;
  //
  vector <Double_t> ztar_foil;
  Int_t ndelcut;
  vector<Double_t > delcut;
  vector<Double_t > delwidth;
  if (file_optics.is_open()) {
    //
    cout << " Open file = " << OpticsFile << endl;
    while (RunNum!=nrun  ) {
      temp.ReadToDelim(file_optics,',');
      cout << temp << endl;
      if (temp.Atoi() == nrun) {
	RunNum = temp.Atoi();
      } else {
	temp.ReadLine(file_optics);
      }
    }
    if (RunNum==nrun) {
      temp.ReadToDelim(file_optics,',');
      OpticsID = temp;
      temp.ReadToDelim(file_optics,',');
      CentAngle = temp.Atof();
      temp.ReadToDelim(file_optics,',');
      NumFoil = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      SieveFlag = temp.Atoi();
      temp.ReadToDelim(file_optics);
      ndelcut = temp.Atoi();
      for (Int_t nf=0;nf<NumFoil-1;nf++) {
        temp.ReadToDelim(file_optics,',');
	ztar_foil.push_back(temp.Atof());
      }
      temp.ReadToDelim(file_optics);
      ztar_foil.push_back(temp.Atof());
      for (Int_t nd=0;nd<ndelcut-1;nd++) {
        temp.ReadToDelim(file_optics,',');
	delcut.push_back(temp.Atof());
      }
      temp.ReadToDelim(file_optics);
      delcut.push_back(temp.Atof());
      for (Int_t nw=0;nw<ndelcut-1;nw++) {
	temp.ReadToDelim(file_optics,',');
	delwidth.push_back(temp.Atof());
      }
      temp.ReadToDelim(file_optics);
      delwidth.push_back(temp.Atof());
    }
  } else {
    cout << " No file = " << OpticsFile << endl;    
  }
  cout << RunNum << " " << OpticsID << " " << CentAngle << " " << NumFoil << " " << SieveFlag << endl;
  if (NumFoil==0) return;
  //
  TString inputroot;
  TFile *fhistroot;
  TString outputpdf;
  inputroot=Form("hist/Optics_%s_%d_hist.root",OpticsID.Data(),FileID);
  //inputroot=Form("hist/Optics_sim_gen_4_-1_hist.root");

  outputpdf = Form("plots/Optics_%s_%d_ytar_cuts",OpticsID.Data(),FileID);
  cout << " infile root = " << inputroot << endl;
  fhistroot =  new TFile(inputroot);
  //
  TH2F *fhist;
  fhist = (TH2F*)fhistroot->Get("hYtarDelta");
  TString YtarDeltaCutFile;
  YtarDeltaCutFile=Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
  TFile *fYtarDeltaCut;
  vector <TCutG*> ytar_delta_cut;
  fYtarDeltaCut = new TFile(YtarDeltaCutFile);
  cout << "Ytar Cut file = " << YtarDeltaCutFile << endl;
  for (Int_t nc=0;nc<NumFoil;nc++) {
    fYtarDeltaCut->cd();
    TCutG* tempcut = (TCutG*)gROOT->FindObject(Form("delta_vs_ytar_cut_foil%d",nc));
    if (tempcut) {
      Int_t npt = tempcut->GetN();
      cout << "hYtarDelta_cut = " << nc << " npts = " << npt << endl;
      ytar_delta_cut.push_back(tempcut);
    } else {
      cout << " No hYtarDelta_cut = " << nc << endl;
    }
  }
  //
  vector<vector<TH2F*> > hYsXs_DelCut;
  vector<TH2F*> temp2d;
  TH2F* th;
  for  (Int_t nc=0;nc<NumFoil;nc++) {
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      th =  (TH2F*)fhistroot->Get(Form("hYsXs_Foil_%d_DelCut_%d",nc,nd));
      if (!th) {
	cout << " no hist : "<< Form("hYsXs_Foil_%d_DelCut_%d",nc,nd) << endl;
	return;
      }
      temp2d.push_back(th);
    }	  
    hYsXs_DelCut.push_back(temp2d);
    temp2d.clear();	
  }
  //
  TLine* xs_line[11];
  TText* xs_text[11];
  TLine* ys_line[11];
  TText* ys_text[11];
  for (Int_t nys=0;nys<11;nys++) {
    Double_t pos=nys*2.5-2.5*5;
    Double_t ypos=nys*1.64-1.64*5;
    xs_line[nys]= new TLine(-12.,pos,12.,pos);
    xs_text[nys]= new TText(-14,pos,Form("%d",nys));
    xs_text[nys]->SetTextColor(2);
    xs_line[nys]->SetLineColor(2);
    xs_line[nys]->SetLineWidth(1);
    ys_line[nys]= new TLine(ypos,-15.,ypos,15);
    ys_text[nys]= new TText(ypos,-17,Form("%d",nys));
    ys_text[nys]->SetTextColor(2);
    ys_line[nys]->SetLineColor(2);
    ys_line[nys]->SetLineWidth(1);
  }
  /*
    TCanvas* can2d[NumFoil][ndelcut];
    for  (Int_t nc=0;nc<NumFoil;nc++) {
    for  (Int_t nd=0;nd<ndelcut;nd++) {
    can2d[nc][nd] = new TCanvas(Form("Can2d_%d_%d",nc,nd),Form("Foil %d Del %d",nc,nd), 700,700);
    can2d[nc][nd]->Divide(1,1);
    can2d[nc][nd]->cd(1);
    hYsXs_DelCut[nc][nd]->Draw("colz");
    hYsXs_DelCut[nc][nd]->SetMinimum(1);
    for (Int_t nys=0;nys<11;nys++) { xs_line[nys]->Draw();}
    for (Int_t nys=0;nys<11;nys++) { xs_text[nys]->Draw();}
    for (Int_t nys=0;nys<11;nys++) { ys_line[nys]->Draw();}
    for (Int_t nys=0;nys<11;nys++) { ys_text[nys]->Draw();}
    TString end = ".pdf";
    if (nc==0 && nd==0) end=".pdf(";
    if (nc==NumFoil-1&& nd==ndelcut-1) end=".pdf)";
    can2d[nc][nd]->Print(outputpdf+end);
    }}
  */
  //
  TCanvas *canCut = new TCanvas("canCut","canCut",700,700);
  canCut->Divide(1,1);
  canCut->cd(1);
  fhist->Draw("colz");
  gPad->SetLogz(); //optional   
  
  //fhist->SetMinimum(100);
  TCutG*t;
  for (Int_t nf=0;nf<NumFoil;nf++) {
    if (ytar_delta_cut[nf]) ytar_delta_cut[nf]->Draw("same");
    if (ytar_delta_cut[nf])	ytar_delta_cut[nf]->SetLineColor(1);
    
  } 
  outputpdf = Form("plots/Optics_%s_%d_ytar_cuts",OpticsID.Data(),FileID);
  canCut->Print(outputpdf+".pdf");
  //
}
