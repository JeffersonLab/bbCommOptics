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


void set_ytar_delta_cuts(Int_t nrun=1813,Int_t FileID=-2) {
  
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
  Int_t CentAngle=0.;
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
   inputroot=Form("hist/Optics_%s_%d_hist.root",OpticsID.Data(),FileID);
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);

   TH2F *fhist;
   fhist = (TH2F*)fhistroot->Get("hYtarDelta");

 
  vector <TString> hname_cut;
  for (Int_t nf=0;nf<NumFoil;nf++) {
    TString temp = Form("delta_vs_ytar_cut_foil%d", nf);
    hname_cut.push_back(temp);
  }
 

  
  
  //Set Name of ROOTfile containing polygon Cuts
  TString outCutFile;
  outCutFile=Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
     cout << " cutfile = " << outCutFile << endl;

  //------------------ Flag 1 ----> Set the Polygon Cut ----------------------

  

    TCanvas *histView_Cut; 
    histView_Cut= new TCanvas("histView_Cut","cut",900,500);
    histView_Cut->Divide(1,1);

    //Loop over focal plane variables
	histView_Cut->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogz(); //optional   
	
	fhistroot->cd();
	
	fhist->Draw("colz");
	fhist->SetMinimum(1);

	
    cout << "outCutFile =  " << outCutFile << endl;
    TFile f(outCutFile,"UPDATE");
   
    for (Int_t nf=0;nf<NumFoil;nf++) {

      TCutG*t;

      Int_t nloop=0;
      while( nloop!=-1) {
        t=(TCutG*)gROOT->FindObject(hname_cut[nf]);
	fhistroot->cd();
	fhist->Draw("colz");
        if(!t){cout << "No cut  =  " <<hname_cut[nf]  << endl;}
 	if(t) {
	  cout << " draw cut = " <<hname_cut[nf] << endl;
	  f.cd();
 		 t->Draw("same");
		 t->SetLineColor(1);
	}
	histView_Cut->Update();
	cout <<" Enter 1 to set cut for foil # " << nf << " (or -1 next foil, -10 delete cut) "  << endl;
	cin >> nloop ;
	if(nloop==-10) {
	  f.cd();
	  f.Delete(hname_cut[nf]+";1");
	  t=(TCutG*)gROOT->FindObject(hname_cut[nf]+";1");
	  gROOT->Remove(t);
	  t=(TCutG*)gROOT->FindObject(hname_cut[nf]+";1");
	  if (!t) cout << " delete cut = " <<hname_cut[nf]  << endl;
	  if (t) cout << " delete cut? = " <<hname_cut[nf]  << endl;
	  f.Write("",TObject::kOverwrite);
	  cout << nf << " delete cut = " <<hname_cut[nf]  << endl;
	}
	
	if (nloop!=-1 && nloop!=-10) {
	  cout << " set cut foil = " << nf << endl;
	  f.cd();
	  f.Delete(hname_cut[nf]+";1");
	  f.Write("",TObject::kOverwrite);
	  t=(TCutG*)gROOT->FindObject(hname_cut[nf]+";1");
	  gROOT->Remove(t);
	  t=(TCutG*)gROOT->FindObject(hname_cut[nf]);
	  if (!t) cout << " delete cut = " <<hname_cut[nf]  << endl;
	  if (t) cout << " delete cut? = " <<hname_cut[nf]  << endl;
	  fhistroot->cd();
	  TCutG*cutg=(TCutG*)gPad->WaitPrimitive("CUTG","CutG");
	  histView_Cut->Update();
	  TCutG*tmpg= (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
	  TCutG*mycutg=(TCutG*)(tmpg->Clone(hname_cut[nf]));
	  f.cd();
	  mycutg->Write("",TObject::kOverwrite);
	  mycutg->Print();
	  mycutg->Draw();
	  histView_Cut->Update();
	}
      } // while
      
      
    }

}


