#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

void plot_ytar_diff(Int_t nrun=1814,Int_t FileID=-2){
   const int nxsieve = 13;
  const int nysieve = 7;
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
  const Int_t Number=3;
  Double_t Red[Number] = { 1.0,0.0,0.0};
  Double_t Blue[Number] = { 1.0,0.0,1.0};
  Double_t Green[Number] = { 0.0,1.0,0.0};
 Double_t Len[Number] = { 0.0,.5,1.0};
 Int_t nb=50;
 TColor::CreateGradientColorTable(Number,Len,Red,Green,Blue,nb);
  //
 TObjArray HList(0);
 //  
 TString  tnrun=Form("%d",nrun);
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
   TString outputhist;
   TString outputpdf;
   TString outputpdf1;
   TString outputpdf2;
   inputroot= Form("hist/Optics_%s_%d_fit_tree.root",OpticsID.Data(),FileID);
   outputpdf= Form("plots/Optics_%s_%d_ydiff",OpticsID.Data(),FileID);
   outputpdf1= Form("plots/Optics_%s_%d_ydiff_xs",OpticsID.Data(),FileID);
   outputpdf2= Form("plots/Optics_%s_%d_ydiff_delta",OpticsID.Data(),FileID);
   outputhist= Form("hist/Optics_%s_%d_ydiff_hist.root",OpticsID.Data(),FileID);
    TFile *hroot = new TFile(inputroot);;
   TTree *FitTree = (TTree*) hroot->Get("TFit");
   Double_t  ys,xtar,ztar,xptar,yptar,ytar,delta,xptarT,yptarT,ytarT,ztarT;
   Double_t ysieveT,ysieve,xsieve,xsieveT,pinvtheta;
   FitTree->SetBranchAddress("ys",&ysieve);
   FitTree->SetBranchAddress("xs",&xsieve);
   FitTree->SetBranchAddress("ysT",&ysieveT);
   FitTree->SetBranchAddress("xsT",&xsieveT);
   FitTree->SetBranchAddress("xtar",&xtar);
   FitTree->SetBranchAddress("ztar",&ztar);
   FitTree->SetBranchAddress("xptar",&xptar);
   FitTree->SetBranchAddress("yptar",&yptar);
   FitTree->SetBranchAddress("ytar",&ytar);
   FitTree->SetBranchAddress("xptarT",&xptarT);
   FitTree->SetBranchAddress("yptarT",&yptarT);
   FitTree->SetBranchAddress("ytarT",&ytarT);
   FitTree->SetBranchAddress("ztarT",&ztarT);
   vector<double> *xpfp = 0;
   FitTree->SetBranchAddress("xpfp",&xpfp);
   vector<double> *ypfp = 0;
   FitTree->SetBranchAddress("ypfp",&ypfp);
   vector<double> *xfp = 0;
   FitTree->SetBranchAddress("xfp",&xfp);
   vector<double> *yfp = 0;
   FitTree->SetBranchAddress("yfp",&yfp);
   FitTree->SetBranchAddress("pinvtheta",&pinvtheta);
   
   //
   
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
   
   vector <Double_t> ys_cent;
   for (Int_t nys=0;nys<7;nys++) {
     Double_t pos=nys*0.0381-0.0381*3;
     ys_cent.push_back(pos);
   }

  
   //
   //
   vector<vector<vector<TH1F*> > > hYDiff;
   vector<vector<TH1F*> > hYDiff_ztar;
   vector<vector<TH2F*> > hYDiff_YpTrue;
   vector<vector<TGraphErrors*> > gYDiff_YpTrue;
   vector<TGraphErrors*>  gYDiff_ztar_delta;
   vector<TH1F*> hYtar;
   vector<TH1F*> hZtar;
   vector<vector<vector<Double_t > > > hYDiff_mean;
   vector<vector<vector<Double_t > > > hYDiff_sigma;
   vector<vector<vector<Double_t > > > hYDiff_cent;
   vector<vector<vector<Double_t > > > hYDiff_centerr;
   vector<vector<Double_t > >  hYDiff_ztar_mean;
   vector<vector<Double_t > >  hYDiff_ztar_sigma;
   vector<vector<Double_t > >  hYDiff_ztar_cent;
   vector<vector<Double_t > >  hYDiff_ztar_centerr;
   
   hYDiff_YpTrue.resize(NumFoil);
   gYDiff_ztar_delta.resize(NumFoil);
   gYDiff_YpTrue.resize(NumFoil);
   hYDiff_ztar.resize(NumFoil);
   hYDiff.resize(NumFoil);
   hYDiff_mean.resize(NumFoil);
   hYDiff_sigma.resize(NumFoil);
   hYDiff_cent.resize(NumFoil);
   hYDiff_centerr.resize(NumFoil);
   hYDiff_ztar_mean.resize(NumFoil);
   hYDiff_ztar_sigma.resize(NumFoil);
   hYDiff_ztar_cent.resize(NumFoil);
   hYDiff_ztar_centerr.resize(NumFoil);
   hYtar.resize(NumFoil);
   hZtar.resize(NumFoil);
   for (Int_t nf=0;nf<NumFoil;nf++) {
     hYDiff_ztar[nf].resize(ndelcut);
     hYDiff_YpTrue[nf].resize(ndelcut);
     gYDiff_YpTrue[nf].resize(ndelcut);
     hYDiff[nf].resize(ndelcut);
     hYDiff_mean[nf].resize(ndelcut);
     hYDiff_sigma[nf].resize(ndelcut);
     hYDiff_cent[nf].resize(ndelcut);
     hYDiff_centerr[nf].resize(ndelcut);
     hYDiff_ztar_mean[nf].resize(ndelcut);
     hYDiff_ztar_sigma[nf].resize(ndelcut);
     hYDiff_ztar_cent[nf].resize(ndelcut);
     hYDiff_ztar_centerr[nf].resize(ndelcut);
     for (Int_t nd=0;nd<ndelcut;nd++) {
       hYDiff[nf][nd].resize(nysieve);
       hYDiff_mean[nf][nd].resize(nysieve);
       hYDiff_sigma[nf][nd].resize(nysieve);
       hYDiff_cent[nf][nd].resize(nysieve);
       hYDiff_centerr[nf][nd].resize(nysieve);
     }}
   //  
   for (Int_t nf=0;nf<NumFoil;nf++) {
     hYtar[nf] = new TH1F(Form("hYtar_%d",nf),Form("Run %s Ztar %3.2f ; Ytar",tnrun.Data(),ztar_foil[nf]),100,-0.3,0.3);
     HList.Add(hYtar[nf]);
     hZtar[nf] = new TH1F(Form("hZtar_%d",nf),Form("Run %s Ztar %3.2f ; Ztar",tnrun.Data(),ztar_foil[nf]),100,-0.3,0.3);
     HList.Add(hZtar[nf]);
     for (Int_t nd=0;nd<ndelcut;nd++) {
       Double_t DelCent=delcut[nd];//(delcut[nd+1]+delcut[nd])/2;
       hYDiff_YpTrue[nf][nd]  = new TH2F(Form("hYDiff_YTrue_%d_DelCut_%d",nf,nd),Form("Run %s Ztar %3.2f Cut %3.1f; Yptar_true (rad); Ytar - Ytar_true (mr) ",tnrun.Data(),ztar_foil[nf],DelCent),90,-.2,.2,50,-0.2,0.2);
       HList.Add(hYDiff_YpTrue[nf][nd]);
       hYDiff_ztar[nf][nd]  = new TH1F(Form("hYDiff_%d_DelCut_%d",nf,nd),Form("Run %s Ztar %3.2f Cut %3.1f; Ytar - Ytar_true (cm)",tnrun.Data(),ztar_foil[nf],DelCent),50,-0.2,0.2);
       HList.Add(hYDiff_ztar[nf][nd]);
       for (Int_t ns=0;ns<nysieve;ns++) {
	 hYDiff[nf][nd][ns]  = new TH1F(Form("hYDiff_%d_DelCut_%d_Xs_%d",nf,nd,ns),Form("Run %s Ztar %3.2f Ys = %d Cut %3.1f; Ytar - Ytar_true (cm)",tnrun.Data(),ztar_foil[nf],ns,DelCent),50,-0.3,0.3);
	 HList.Add(hYDiff[nf][nd][ns]);
       }  }}
   //
    
 //
   Long64_t nentries = FitTree->GetEntries();
   cout << " start loop " << nentries << endl;
   for (int i = 0; i < nentries; i++) {
     FitTree->GetEntry(i);
     for (Int_t nf=0;nf<NumFoil;nf++) {
       if (abs(ztarT-ztar_foil[nf])<0.02) { 
	 for (Int_t nd=0;nd<ndelcut;nd++) {
	   if ( pinvtheta >=delcut[nd]-delwidth[nd] && pinvtheta <delcut[nd]+delwidth[nd]) {
	     hYtar[nf]->Fill(ytar);
	     hZtar[nf]->Fill(ztar);
	     hYDiff_YpTrue[nf][nd]->Fill(yptarT,(ytar-ytarT));
	     hYDiff_ztar[nf][nd]->Fill((ytar-ytarT));
	     for (Int_t ns=0;ns<nysieve;ns++) {
	       if ( abs(ysieveT-ys_cent[ns])<0.02) hYDiff[nf][nd][ns]->Fill((ytar-ytarT));
	     }
	   }
	 }
       }
     }		
     //	
   }
   //
   TCanvas* can2d[NumFoil];
   Int_t colind=1;
   for  (Int_t nc=0;nc<NumFoil;nc++) {
     can2d[nc] = new TCanvas(Form("Can2d_%d",nc),Form("Foil %d",nc), 700,700);
     //can2d[nc]->Divide(2,5);
     for  (Int_t nd=0;nd<ndelcut;nd++) {
       //can2d[nc]->cd(nd+1);
       //	  hYDiff_YpTrue[nc][nd]->Draw("colz");
       Double_t DelCent=delcut[nd];//(delcut[nd+1]+delcut[nd])/2;
       hYDiff_ztar[nc][nd]->Draw("colz");
       hYDiff_ztar_mean[nc][nd]= hYDiff_ztar[nc][nd]->GetMean();
       hYDiff_ztar_sigma[nc][nd]= hYDiff_ztar[nc][nd]->GetRMS();
       hYDiff_ztar_cent[nc][nd]= DelCent;
       hYDiff_ztar_centerr[nc][nd]= 0.0001;
     }
     gYDiff_ztar_delta[nc] = new TGraphErrors(ndelcut,&hYDiff_ztar_cent[nc][0],&hYDiff_ztar_mean[nc][0],&hYDiff_ztar_centerr[nc][0],&hYDiff_ztar_sigma[nc][0]);
     if (colind==5) colind++;
     gYDiff_ztar_delta[nc]->SetMarkerColor(colind++);
     gYDiff_ztar_delta[nc]->SetMarkerStyle(nc+20);
     TString end = ".pdf";
     if (nc==0) end=".pdf(";
     if (nc==NumFoil-1) end=".pdf)";
     can2d[nc]->Print(outputpdf+end);
     
   }  //
   //
   TCanvas* candelta;
   TMultiGraph *mgrdelta;
   TLegend *legdelta;
   candelta = new TCanvas("candelta","candelta", 700,700);
   candelta->Divide(1,1);
   candelta->cd(1);
   legdelta = new TLegend(.79,.65,.99,.95,"");
   mgrdelta=new TMultiGraph();  
   for  (Int_t nf=0;nf<NumFoil;nf++) {
     mgrdelta->Add(gYDiff_ztar_delta[nf]);
     legdelta->AddEntry(gYDiff_ztar_delta[nf],Form("Ztar = %3.1f",ztar_foil[nf]),"p");
   }
   mgrdelta->SetTitle(Form("SHMS Angle = %4.2f; Delta ; Ytar -Y_true (cm)",CentAngle));
   mgrdelta->SetMinimum(-0.2);
   mgrdelta->SetMaximum(+0.2);
   mgrdelta->Draw("AP");
   legdelta->Draw();
   TString end = ".pdf";
   candelta->Print(outputpdf2+end);
   
   //	TCanvas* candel[NumFoil][ndelcut];
   Double_t Y_offset=0.; // Offset to get delta = 0 ztar=0 to about zero.
   for  (Int_t nf=0;nf<NumFoil;nf++) {
     Int_t colind = 0;
     for  (Int_t nd=0;nd<ndelcut;nd++) {
       //	  candel[nf][nd] = new TCanvas(Form("Candel_%d_%d",nf,nd),Form("Candel_%d_%d",nf,nd), 700,700);
       //	  candel[nf][nd]->Divide(2,5);
       for (Int_t ns=0;ns<7;ns++) {
	 //   candel[nf][nd]->cd(ns+1);
	 //  hYDiff[nf][nd][ns]->Draw();
	 hYDiff_mean[nf][nd][ns]= hYDiff[nf][nd][ns]->GetMean()-Y_offset;
	 hYDiff_sigma[nf][nd][ns]= hYDiff[nf][nd][ns]->GetRMS();
	 hYDiff_cent[nf][nd][ns]= ys_cent[ns];
	 hYDiff_centerr[nf][nd][ns]= 0.0001;
       }
       gYDiff_YpTrue[nf][nd] = new TGraphErrors(7,&hYDiff_cent[nf][nd][0],&hYDiff_mean[nf][nd][0],&hYDiff_centerr[nf][nd][0],&hYDiff_sigma[nf][nd][0]);
       colind++;
       if (colind==5) colind++;
       if (colind==8) colind=1;
       gYDiff_YpTrue[nf][nd]->SetMarkerColor(colind);
       gYDiff_YpTrue[nf][nd]->SetMarkerStyle(nd+20);
     }}	  
   gStyle->SetPadRightMargin(0.2);
   
   //
   TCanvas* candel[NumFoil];
   TMultiGraph *mgr[NumFoil];
   TLegend *leg[NumFoil];
   for  (Int_t nf=0;nf<NumFoil;nf++) {
     candel[nf] = new TCanvas(Form("Candel_%d",nf),Form("Candel_%d",nf), 700,700);
     candel[nf]->Divide(1,1);
     leg[nf] = new TLegend(.79,.65,.99,.95,"");
     mgr[nf]=new TMultiGraph();  
     for  (Int_t nd=0;nd<ndelcut;nd++) {
       mgr[nf]->Add(gYDiff_YpTrue[nf][nd]);
       Double_t DelCent=delcut[nd];//(delcut[nd+1]+delcut[nd])/2;
       leg[nf]->AddEntry(gYDiff_YpTrue[nf][nd],Form("Cut = %3.2f",DelCent),"p");
     }
     candel[nf]->cd(1);
     mgr[nf]->SetTitle(Form("Ztar = %4.2f BB Angle = %4.2f; Y_sieve (cm); Ytar -Y_true (cm)",ztar_foil[nf],CentAngle));
     mgr[nf]->SetMinimum(-0.2);
     mgr[nf]->SetMaximum(+0.2);
     mgr[nf]->Draw("AP");
     leg[nf]->Draw();
     TString end = ".pdf";
     if (nf==0) end=".pdf(";
     if (nf==NumFoil-1) end=".pdf)";
     candel[nf]->Print(outputpdf1+end);
   }
   //
   TFile hsimc(outputhist,"recreate");
   HList.Write();
   
//
}
