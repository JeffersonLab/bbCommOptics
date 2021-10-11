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

void plot_xptar_diff(Int_t nrun=1814,Int_t FileID=-2){
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
  const int nxsieve = 13;
  const int nysieve = 7;
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
   inputroot= Form("hist/Optics_%s_%d_fit_tree.root",OpticsID.Data(),FileID);
   outputpdf= Form("plots/Optics_%s_%d_xpdiff",OpticsID.Data(),FileID);
   outputpdf1= Form("plots/Optics_%s_%d_xpdiff_xs",OpticsID.Data(),FileID);
   outputhist= Form("hist/Optics_%s_%d_xpdiff_hist.root",OpticsID.Data(),FileID);
    TFile *hroot = new TFile(inputroot);;
   TTree *FitTree = (TTree*) hroot->Get("TFit");
   Double_t  ys,xtar,ztar,xptar,yptar,ytar,delta,xptarT,yptarT,ytarT,ztarT;
   Double_t ysieveT,ysieve, xsieve, xsieveT,pinvtheta;
   FitTree->SetBranchAddress("ys",&ysieve);
   FitTree->SetBranchAddress("ysT",&ysieveT);
    FitTree->SetBranchAddress("xs",&xsieve);
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

 vector <Double_t> ys_cent;
  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=nys*0.0381-0.0381*3;//old sieve
    ys_cent.push_back(pos);
  }

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
 
 /*
   vector <Double_t> ys_cent;
   for (Int_t nys=0;nys<7;nys++) {
     Double_t pos=nys*0.0381-0.0381*3;
     ys_cent.push_back(pos);
   }
 */
   //
 
 //
 //
  vector<vector<vector<TH1F*> > > hXpDiff;
  vector<vector<vector<vector<TH1F*> > > >hXpDiff_Ys;
  vector<vector<TH2F*> > hXpDiff_XpTrue;
  vector<vector<TGraphErrors*> > gXpDiff_XpTrue;
  vector<vector<vector<TGraphErrors*> > > gXpDiff_Ys_XpTrue;
  vector<TH1F*> hYtar;
  vector<TH1F*> hZtar;
  vector<vector<vector<Double_t > > > hXpDiff_mean;
  vector<vector<vector<Double_t > > > hXpDiff_sigma;
  vector<vector<vector<Double_t > > > hXpDiff_cent;
  vector<vector<vector<Double_t > > > hXpDiff_centerr;
  vector<vector<vector<vector<Double_t > > > > hXpDiff_Ys_mean;
  vector<vector<vector<vector<Double_t > > > > hXpDiff_Ys_sigma;
  vector<vector<vector<vector<Double_t > > > > hXpDiff_Ys_cent;
  vector<vector<vector<vector<Double_t > > > > hXpDiff_Ys_centerr;
  
  hXpDiff_XpTrue.resize(NumFoil);
  gXpDiff_XpTrue.resize(NumFoil);
  gXpDiff_Ys_XpTrue.resize(NumFoil);
  hXpDiff_Ys.resize(NumFoil);
  hXpDiff.resize(NumFoil);
  hXpDiff_mean.resize(NumFoil);
  hXpDiff_sigma.resize(NumFoil);
  hXpDiff_cent.resize(NumFoil);
  hXpDiff_centerr.resize(NumFoil);
  hXpDiff_Ys_mean.resize(NumFoil);
  hXpDiff_Ys_sigma.resize(NumFoil);
  hXpDiff_Ys_cent.resize(NumFoil);
  hXpDiff_Ys_centerr.resize(NumFoil);
  hYtar.resize(NumFoil);
  hZtar.resize(NumFoil);
   for (Int_t nf=0;nf<NumFoil;nf++) {
  hXpDiff_XpTrue[nf].resize(ndelcut);
  gXpDiff_XpTrue[nf].resize(ndelcut);
  gXpDiff_Ys_XpTrue[nf].resize(ndelcut);
  hXpDiff[nf].resize(ndelcut);
  hXpDiff_Ys[nf].resize(ndelcut);
  hXpDiff_mean[nf].resize(ndelcut);
  hXpDiff_sigma[nf].resize(ndelcut);
  hXpDiff_cent[nf].resize(ndelcut);
  hXpDiff_centerr[nf].resize(ndelcut);
  hXpDiff_Ys_mean[nf].resize(ndelcut);
  hXpDiff_Ys_sigma[nf].resize(ndelcut);
  hXpDiff_Ys_cent[nf].resize(ndelcut);
  hXpDiff_Ys_centerr[nf].resize(ndelcut);
   for (Int_t nd=0;nd<ndelcut;nd++) {
      hXpDiff[nf][nd].resize(nxsieve);
      hXpDiff_Ys[nf][nd].resize(nxsieve);
      hXpDiff_mean[nf][nd].resize(nxsieve);
      hXpDiff_sigma[nf][nd].resize(nxsieve);
      hXpDiff_cent[nf][nd].resize(nxsieve);
      hXpDiff_centerr[nf][nd].resize(nxsieve);
      hXpDiff_Ys_mean[nf][nd].resize(nxsieve);
      hXpDiff_Ys_sigma[nf][nd].resize(nxsieve);
      hXpDiff_Ys_cent[nf][nd].resize(nxsieve);
      hXpDiff_Ys_centerr[nf][nd].resize(nxsieve);
      gXpDiff_Ys_XpTrue[nf][nd].resize(nxsieve);
      for (Int_t ny=0;ny<nxsieve;ny++) {
	hXpDiff_Ys[nf][nd][ny].resize(nxsieve);
      hXpDiff_Ys_mean[nf][nd][ny].resize(nxsieve);
      hXpDiff_Ys_sigma[nf][nd][ny].resize(nxsieve);
      hXpDiff_Ys_cent[nf][nd][ny].resize(nxsieve);
      hXpDiff_Ys_centerr[nf][nd][ny].resize(nxsieve);
      }

   }}
   //  
    for (Int_t nf=0;nf<NumFoil;nf++) {
      hYtar[nf] = new TH1F(Form("hYtar_%d",nf),Form("Run %s Ztar %3.2f ; Ytar",tnrun.Data(),ztar_foil[nf]),100,-0.3,0.3);
      HList.Add(hYtar[nf]);
      hZtar[nf] = new TH1F(Form("hZtar_%d",nf),Form("Run %s Ztar %3.2f ; Ztar",tnrun.Data(),ztar_foil[nf]),100,-0.3,0.3);
      HList.Add(hZtar[nf]);
      for (Int_t nd=0;nd<ndelcut;nd++) {
	Double_t DelCent=delcut[nd];//(delcut[nd+1]+delcut[nd])/2;
	hXpDiff_XpTrue[nf][nd]  = new TH2F(Form("hXpDiff_XpTrue_%d_DelCut_%d",nf,nd),Form("Run %s Ztar %3.2f Cut %3.1f; Xptar_true (rad); Xptar_true - Xptar (r) ",tnrun.Data(),ztar_foil[nf],DelCent),90,-0.3,0.3,80,-0.3,0.3);
        HList.Add(hXpDiff_XpTrue[nf][nd]);
	for (Int_t nxs=0;nxs<nxsieve;nxs++) {
	  hXpDiff[nf][nd][nxs]  = new TH1F(Form("hXpDiff_%d_DelCut_%d_Xs_%d",nf,nd,nxs),Form("Run %s Ztar %3.2f Xs = %d Cut %3.1f; Xptar_true - Xptar (r)",tnrun.Data(),ztar_foil[nf],nxs,DelCent),80,-0.3,0.3);
	  HList.Add(hXpDiff[nf][nd][nxs]);
	  for (Int_t nys=0;nys<nxsieve;nys++) {
	    hXpDiff_Ys[nf][nd][nxs][nys]  = new TH1F(Form("hXpDiff_%d_DelCut_%d_Xs_%d_Ys_%d",nf,nd,nxs,nys),Form("Run %s Ztar %3.3f Xs = %d  Ys = %d Cut %3.1f; Xptar_true - Xptar (r)",tnrun.Data(),ztar_foil[nf],nxs,nys,DelCent),80,-0.3,0.3);
	    HList.Add(hXpDiff_Ys[nf][nd][nxs][nys]);
    }
    }	
    }  
    }
    //
    
 //
Long64_t nentries = FitTree->GetEntries();
 cout << " start loop " << nentries << endl;
	for (int i = 0; i < nentries; i++) {
      		FitTree->GetEntry(i);
    for (Int_t nf=0;nf<NumFoil;nf++) {
      if (abs(ztarT-ztar_foil[nf])<0.02) { 
       for (Int_t nd=0;nd<ndelcut;nd++) {
	 if ( delta >=delcut[nd]-delwidth[nd] && delta <delcut[nd]+delwidth[nd]) {
	hYtar[nf]->Fill(ytar);
	hZtar[nf]->Fill(ztar);
	   hXpDiff_XpTrue[nf][nd]->Fill(xptarT,(xptar-xptarT));
           for (Int_t nxs=0;nxs<nxsieve;nxs++) {
	     if ( abs(xsieveT-xs_cent[nxs])<0.02) {
	       hXpDiff[nf][nd][nxs]->Fill((xptar-xptarT));
               for (Int_t nys=0;nys<nysieve;nys++) {
		 if ( abs(ysieveT-ys_cent[nys])<.02)hXpDiff_Ys[nf][nd][nxs][nys]->Fill((xptar-xptarT));
	       }
	     }
	   }
	 }
       }
      }
    }		
    //	
	}
	/*
	  TCanvas* can2d[NumFoil];
	  for  (Int_t nc=0;nc<NumFoil;nc++) {
	  can2d[nc] = new TCanvas(Form("Can2d_%d",nc),Form("Foil %d",nc), 700,700);
	  can2d[nc]->Divide(2,5);
	  for  (Int_t nd=0;nd<ndelcut;nd++) {
	  can2d[nc]->cd(nd+1);
	  hXpDiff_XpTrue[nc][nd]->Draw("colz");
	  }
	  TString end = ".pdf";
	  if (nc==0) end=".pdf(";
	  if (nc==NumFoil-1) end=".pdf)";
	  can2d[nc]->Print(outputpdf+end);
	  
	  }  //
	*/
	//	TCanvas* candel[NumFoil][ndelcut];
	Double_t Xp_offset=0.; // Offset to get delta = 0 ztar=0 to about zero.
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	  Int_t colind = 0;
	  for  (Int_t nd=0;nd<ndelcut;nd++) {
	    //	  candel[nf][nd] = new TCanvas(Form("Candel_%d_%d",nf,nd),Form("Candel_%d_%d",nf,nd), 700,700);
	    //	  candel[nf][nd]->Divide(2,5);
	    for (Int_t nxs=0;nxs<nxsieve;nxs++) {
	      //   candel[nf][nd]->cd(nxs+1);
	      //  hXpDiff[nf][nd][nxs]->Draw();
	      hXpDiff_mean[nf][nd][nxs]= hXpDiff[nf][nd][nxs]->GetMean()-Xp_offset;
	      hXpDiff_sigma[nf][nd][nxs]= hXpDiff[nf][nd][nxs]->GetRMS();
	      hXpDiff_cent[nf][nd][nxs]= xs_cent[nxs];
	      hXpDiff_centerr[nf][nd][nxs]= 0.0001;
	      
	      gXpDiff_XpTrue[nf][nd] = new TGraphErrors(nxsieve,&hXpDiff_cent[nf][nd][0],&hXpDiff_mean[nf][nd][0],&hXpDiff_centerr[nf][nd][0],&hXpDiff_sigma[nf][nd][0]);
	      colind++;
	      if (colind==5) colind++;
	      if (colind==8) colind=1;
	      gXpDiff_XpTrue[nf][nd]->SetMarkerColor(colind);
	      gXpDiff_XpTrue[nf][nd]->SetMarkerStyle(nd+20);
	    }
	    //
	    for (Int_t nxs=0;nxs<nxsieve;nxs++) {
	      for (Int_t nys=0;nys<nysieve;nys++) {
		hXpDiff_Ys_mean[nf][nd][nxs][nys]= hXpDiff_Ys[nf][nd][nxs][nys]->GetMean()-Xp_offset;
		hXpDiff_Ys_sigma[nf][nd][nxs][nys]=hXpDiff_Ys[nf][nd][nxs][nys]->GetRMS();
		hXpDiff_Ys_cent[nf][nd][nxs][nys]=ys_cent[nys];
		hXpDiff_Ys_centerr[nf][nd][nxs][nys]=0.0001;
	      }
	      
	    }
	    Int_t colind2=0;
	    for (Int_t nxs=0;nxs<nxsieve;nxs++) {
	      gXpDiff_Ys_XpTrue[nf][nd][nxs] = new TGraphErrors(hXpDiff_Ys_cent[nf][nd][nxs].size(),&hXpDiff_Ys_cent[nf][nd][nxs][0],&hXpDiff_Ys_mean[nf][nd][nxs][0],&hXpDiff_Ys_centerr[nf][nd][nxs][0],&hXpDiff_Ys_sigma[nf][nd][nxs][0]);
	      colind2++;
	      if (colind2==5) colind2++;
	      if (colind2==8) colind2=1;
	      gXpDiff_Ys_XpTrue[nf][nd][nxs]->SetMarkerColor(colind2);
	      gXpDiff_Ys_XpTrue[nf][nd][nxs]->SetMarkerStyle(nxs+20);
	    }
	  }}	  
	gStyle->SetPadRightMargin(0.2);
        
	//
	TCanvas* canFoilDel[NumFoil][ndelcut];
        TMultiGraph *mgrFoilDel[NumFoil][ndelcut];
        TLegend *legFoilDel[NumFoil][ndelcut];
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	  for  (Int_t nd=0;nd<ndelcut;nd++) {
	    canFoilDel[nf][nd] = new TCanvas(Form("CanFoilDe1_%d_%d",nf,nd),Form("CanFoilDe1_%d_%d",nf,nd), 700,700);
	    canFoilDel[nf][nd]->Divide(1,1);
	    legFoilDel[nf][nd] = new TLegend(.79,.65,.99,.95,"");
	    mgrFoilDel[nf][nd]=new TMultiGraph();  
	    for (Int_t nxs=0;nxs<11;nxs++) {
	      mgrFoilDel[nf][nd]->Add(gXpDiff_Ys_XpTrue[nf][nd][nxs]);
	      legFoilDel[nf][nd]->AddEntry(gXpDiff_Ys_XpTrue[nf][nd][nxs],Form("Xs = %3.3f",xs_cent[nxs]),"p");
	    }
	    Double_t DelCent=delcut[nd];//(delcut[nd+1]+delcut[nd])/2;
	    canFoilDel[nf][nd]->cd(1);
	    mgrFoilDel[nf][nd]->SetTitle(Form("Ztar = %4.2f Cut = %3.1f BB Angle = %4.2f; Y_sieve (cm); Xptar -Xp_true (mr)",ztar_foil[nf],DelCent,CentAngle));
	    mgrFoilDel[nf][nd]->SetMinimum(-0.3);
	    mgrFoilDel[nf][nd]->SetMaximum(+0.3);
	    mgrFoilDel[nf][nd]->Draw("AP");
	    legFoilDel[nf][nd]->Draw();
	    TString end = ".pdf";
	    if (nf==0 && nd==0) end=".pdf(";
	    if (nf==NumFoil-1 && nd==ndelcut-1) end=".pdf)";
	    canFoilDel[nf][nd]->Print(outputpdf+end);
	  }}
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
	    mgr[nf]->Add(gXpDiff_XpTrue[nf][nd]);
	    Double_t DelCent=delcut[nd];//(delcut[nd+1]+delcut[nd])/2;
	    leg[nf]->AddEntry(gXpDiff_XpTrue[nf][nd],Form("Cut = %3.1f",DelCent),"p");
	  }
	  candel[nf]->cd(1);
	  mgr[nf]->SetTitle(Form("Ztar = %4.2f BB Angle = %4.2f; X_sieve (cm); Xptar -Xp_true (r)",ztar_foil[nf],CentAngle));
	  mgr[nf]->SetMinimum(-0.3);
	  mgr[nf]->SetMaximum(+0.3);
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
