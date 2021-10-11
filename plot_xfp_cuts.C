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


void plot_xfp_cuts(Int_t nrun=1814,Int_t FileID=-2) {
  const int nxsieve = 13;
  const int nysieve = 7;
  
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
   outputpdf = Form("plots/Optics__%s_%d_xfp_cuts",OpticsID.Data(),FileID);
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
 //
   vector<vector<Double_t> > AxisRange;
   AxisRange.resize(ndelcut);
   TString AxisRangeFileName = "AxisRange_xpfp_xfp.dat";
   ifstream AxisRangeFile(AxisRangeFileName.Data());
   for  (Int_t nd=0;nd<ndelcut;nd++) { AxisRange[nd].resize(4) ;}
  //
  TString temp1;
  if (file_optics.is_open()) {
   for  (Int_t nd=0;nd<ndelcut;nd++) { 
     temp1.ReadToDelim(AxisRangeFile,',');
     AxisRange[nd][0] = temp1.Atof();
     temp1.ReadToDelim(AxisRangeFile,',');
     AxisRange[nd][1] = temp1.Atof();
     temp1.ReadToDelim(AxisRangeFile,',');
     AxisRange[nd][2] = temp1.Atof();
     temp1.ReadToDelim(AxisRangeFile);
     AxisRange[nd][3] = temp1.Atof();
      }
  }
  //
	//
  vector<vector<TH2F*> > hYsXs_DelCut;
  vector<vector<TH2F*> > hXpFpXFp_DelCut;
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
  for  (Int_t nc=0;nc<NumFoil;nc++) {
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      th =  (TH2F*)fhistroot->Get(Form("hXpFpXFp_%d_DelCut_%d",nc,nd));
      if (!th) {
	cout << " no hist : "<< Form("hXpFpXFp_%d_DelCut_%d",nc,nd) << endl;
	return;
      }
      temp2d.push_back(th);
    }	  
    hXpFpXFp_DelCut.push_back(temp2d);
    temp2d.clear();	
  }
  //
  vector<vector<vector<TH2F*> > > hYsXs_DelCut_XpXfpCut;
  hYsXs_DelCut_XpXfpCut.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    hYsXs_DelCut_XpXfpCut[nf].resize(ndelcut);
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      hYsXs_DelCut_XpXfpCut[nf][nd].resize(nxsieve);
    }
  }
  for  (Int_t nc=0;nc<NumFoil;nc++) {
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      for  (Int_t ny=0;ny<nxsieve;ny++) {
	hYsXs_DelCut_XpXfpCut[nc][nd][ny] =  (TH2F*)fhistroot->Get(Form("hYsXs_Foil_%d_DelCut_%d_XFpCut_%d",nc,nd,ny));
      }}}
  //
  TLine* xs_line[nxsieve];
  TText* xs_text[nxsieve];
  
  vector <Double_t> xs_cent{-(0.3+0.0492)+0.0493/cos(18.*3.14/180.),
      -(0.3+0.0492)+(0.0493+0.0492)/cos(18.*3.14/180),
      -(0.3+0.0492)+0.1493/cos(9.*3.14/180.),
      -(0.3+0.0492)+(0.1493+0.0492)/cos(9.*3.14/180.),
      -(0.3+0.0492)+(0.1493+0.0492*2.)/cos(9.*3.14/180.),
      -0.0492,
      0.0,
      0.0492,
      0.3+0.0492-(0.1493+0.0492*2.)/cos(9.*3.14/180.),
      0.3+0.0492-(0.1493+0.0492)/cos(9.*3.14/180.),
      0.3+0.0492-0.1493/cos(9.*3.14/180.),
      0.3+0.0492-(0.0493+0.0492)/cos(18.*3.14/180),
      0.3+0.0492-0.0493/cos(18.*3.14/180.)};
  
  
  for (Int_t nys=0;nys<nxsieve;nys++) {
    Double_t pos=xs_cent[nys];
    xs_line[nys]= new TLine(-0.2,pos,0.2,pos);
    xs_text[nys]= new TText(-0.25,pos,Form("%d",nys));
    xs_text[nys]->SetTextColor(2);
    xs_line[nys]->SetLineColor(2);
    xs_line[nys]->SetLineWidth(1);
  }
  
 
  TLine* ys_line[nysieve];
  TText* ys_text[nysieve];
  
  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=nys*0.0381-0.0381*3;
    ys_line[nys]= new TLine(pos,-0.4,pos,0.4);
    ys_text[nys]= new TText(pos,-0.5,Form("%d",nys));
    ys_text[nys]->SetTextColor(2);
    ys_line[nys]->SetLineColor(2);
    ys_line[nys]->SetLineWidth(1);
  }
  //
	/*
	TCanvas* can2d[NumFoil][ndelcut];
	//	NumFoil=1;
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	    for  (Int_t nd=0;nd<ndelcut;nd++) {
	      can2d[nc][nd] = new TCanvas(Form("Can2d_%d_%d",nc,nd),Form("Foil %d Del %d",nc,nd), 700,700);
	      can2d[nc][nd]->Divide(4,4);
	  can2d[nc][nd]->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogz();    
	  hXpFpXFp_DelCut[nc][nd]->Draw("colz");
	  hXpFpXFp_DelCut[nc][nd]->SetMinimum(1);
	  hXpFpXFp_DelCut[nc][nd]->GetYaxis()->SetRangeUser(AxisRange[nd][0],AxisRange[nd][1]);
	  hXpFpXFp_DelCut[nc][nd]->GetXaxis()->SetRangeUser(AxisRange[nd][2],AxisRange[nd][3]);
	  can2d[nc][nd]->cd(2);
	  hYsXs_DelCut[nc][nd]->Draw("colz");
	  hYsXs_DelCut[nc][nd]->SetMinimum(1);
          for (Int_t nys=0;nys<11;nys++) { xs_line[nys]->Draw();}
          for (Int_t nys=0;nys<11;nys++) { xs_text[nys]->Draw();}
	  Int_t incr=0;
	  	for  (Int_t ny=0;ny<11;ny++) {
		  if (hYsXs_DelCut_XpXfpCut[nc][nd][ny]->Integral()>0) {
		    can2d[nc][nd]->cd(3+incr);
		    hYsXs_DelCut_XpXfpCut[nc][nd][ny]->Draw("colz");
		    hYsXs_DelCut_XpXfpCut[nc][nd][ny]->SetMinimum(1);
		    for (Int_t nys=0;nys<11;nys++) { xs_line[nys]->Draw();}
		    for (Int_t nys=0;nys<11;nys++) { xs_text[nys]->Draw();}
		    incr++;
		  }
		}
	  TString end = ".pdf";
	  if (nc==0 && nd==0) end=".pdf(";
	  if (nc==NumFoil-1 && nd==ndelcut-1) end=".pdf)";
	  can2d[nc][nd]->Print(outputpdf+end);
	    }}
	*/
	//
  outputpdf = Form("plots/Optics__%s_%d_xpfp_xfp_cuts",OpticsID.Data(),FileID);
  TString xpfp_xfp_outCutFile;
  TFile *xpfp_xfp_fcut;
  vector<vector<vector<TCutG*> > > xpfp_xfp_cut;
  vector<vector<vector<Int_t> > > xpfp_xfp_cut_flag;
  xpfp_xfp_cut.resize(NumFoil);
  xpfp_xfp_cut_flag.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    xpfp_xfp_cut[nf].resize(ndelcut);
    xpfp_xfp_cut_flag[nf].resize(ndelcut);
  }
  xpfp_xfp_outCutFile=Form("cuts/XpFpXFp_%s_%d_cut.root",OpticsID.Data(),FileID);
  xpfp_xfp_fcut = new TFile(xpfp_xfp_outCutFile);
  cout << "xpfp_xfp_ Cut file = " << xpfp_xfp_outCutFile << endl;
  xpfp_xfp_fcut->cd();
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      for (Int_t nc=0;nc<nxsieve;nc++) {
	TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	if (tempg)  {
	  //cout << "hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	  xpfp_xfp_cut[nf][nd].push_back(tempg);
	} else {
	  //cout << " No hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
     	  xpfp_xfp_cut[nf][nd].push_back(tempg);
	}
      }}}

  TCanvas* canCut[NumFoil][ndelcut];
  //	NumFoil=1;
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      canCut[nf][nd] = new TCanvas(Form("CanCut_%d_%d",nf,nd),Form("Foil %d Del %d",nf,nd), 700,700);
      canCut[nf][nd]->Divide(2,2);
      canCut[nf][nd]->cd(1);
      gPad->SetGridx();
      gPad->SetGridy();
      //gPad->SetLogz();    
      hXpFpXFp_DelCut[nf][nd]->Draw("colz");
      hXpFpXFp_DelCut[nf][nd]->SetMinimum(1);
      hXpFpXFp_DelCut[nf][nd]->GetYaxis()->SetRangeUser(AxisRange[nd][0],AxisRange[nd][1]);
      hXpFpXFp_DelCut[nf][nd]->GetXaxis()->SetRangeUser(AxisRange[nd][2],AxisRange[nd][3]);
      for (Int_t nc=0;nc<nxsieve;nc++) {
	if (xpfp_xfp_cut[nf][nd][nc]) {
	  xpfp_xfp_cut[nf][nd][nc]->Draw("same");
	  xpfp_xfp_cut[nf][nd][nc]->SetLineColor(2);
	  xpfp_xfp_cut[nf][nd][nc]->SetLineWidth(2);
	  Double_t xcut,ycut;
	  xpfp_xfp_cut[nf][nd][nc]->GetPoint(0,xcut,ycut);
	  TText* ystext = new TText(xcut,ycut,Form("%d",nc));
	  ystext->Draw();
	}
      }

      canCut[nf][nd]->cd(2);
      gPad->SetLogz();    
      hYsXs_DelCut[nf][nd]->Draw("colz");
      hYsXs_DelCut[nf][nd]->SetMinimum(1);
                  	cout<<"no problem yet"<<endl;

      for (Int_t nxs=0;nxs<nxsieve;nxs++) { xs_line[nxs]->Draw();}
      for (Int_t nxs=0;nxs<nxsieve;nxs++) { xs_text[nxs]->Draw();}
      for (Int_t nys=0;nys<nysieve;nys++) { ys_line[nys]->Draw();}
      for (Int_t nys=0;nys<nysieve;nys++) { ys_text[nys]->Draw();}

      canCut[nf][nd]->cd(3);
      gPad->SetLogz();
      /*
      for  (Int_t ny=0;ny<nxsieve;ny++) {
	if (hYsXs_DelCut_XpXfpCut[nf][nd][ny]->Integral()>0) {
	  if (ny==0) {
	    hYsXs_DelCut_XpXfpCut[nf][nd][ny]->Draw("colz");
	    hYsXs_DelCut_XpXfpCut[nf][nd][ny]->SetMinimum(1);
	    hYsXs_DelCut_XpXfpCut[nf][nd][ny]->SetMaximum(300);
	  } else {
	    hYsXs_DelCut_XpXfpCut[nf][nd][ny]->Draw("colz same");
	    hYsXs_DelCut_XpXfpCut[nf][nd][ny]->SetMinimum(1);
	    hYsXs_DelCut_XpXfpCut[nf][nd][ny]->SetMaximum(300);
	  }
	  }
      }
      */

      for (Int_t nys=0;nys<nxsieve;nys++) { xs_line[nys]->Draw("same");}
      for (Int_t nys=0;nys<nxsieve;nys++) { xs_text[nys]->Draw("same");}
      for (Int_t nys=0;nys<nysieve;nys++) { ys_line[nys]->Draw("same");}
      for (Int_t nys=0;nys<nysieve;nys++) { ys_text[nys]->Draw("same");}
      TString end = ".pdf";
      if (nf==0 && nd==0) end=".pdf(";
      if (nf==NumFoil-1 && nd==ndelcut-1) end=".pdf)";
      canCut[nf][nd]->Print(outputpdf+end);
    }}
  //
}
