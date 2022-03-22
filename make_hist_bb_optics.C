#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <TMath.h>
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
#include "TVector3.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

void make_hist_bb_optics(Int_t nrun=1813,Bool_t CutYtarFlag=kTRUE,Bool_t CutYpFpYFpFlag=kTRUE,Bool_t CutXpFpXFpFlag=kTRUE,Int_t FileID=-2){

  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
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
  Int_t ndelcut=-1;
  vector<Double_t > delcut;
  vector<Double_t > delwidth;
  if (file_optics.is_open()) {
    //
    cout << " Open file = " << OpticsFile << endl;
    //while (RunNum!=nrun  ) {
      temp.ReadToDelim(file_optics,',');
      cout << temp << endl;
      if (temp.Atoi() == nrun) {
	RunNum = temp.Atoi();
      } else {
	temp.ReadLine(file_optics);
      }
      //}
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
  CentAngle *= 3.14/180.0;
  
  TString inputroot;
  TString outputhist;
  
  // if (nrun==11107){
    //inputroot=Form("../sim/replayed_simdigtest_2_20211004.root");//shms_replay_matrixopt_%s_%d.root",OpticsID.Data(),FileID);
    //inputroot=Form("Rootfiles/gmn_replayed_11107_merged.root");//stream0_seg0_0.root");  
    //inputroot=Form("Rootfiles/gmn_replayed_11107_11109.root");
    //}
  //else if (nrun==11175){inputroot=Form("Rootfiles/gmn_replayed_11175_11178.root");}
  //else{
  //inputroot=Form("Rootfiles/combined_13674_5foil.root");//}
  inputroot=Form("Rootfiles/combined_13675_4foil_new.root");//}
  //inputroot=Form("files/replay_sbs8_13437_13440.root");//}
  //inputroot=Form("files/gmn_replayed_11966_11985.root");//}
  outputhist=Form("hist/Optics_%s_%d_hist.root",OpticsID.Data(),FileID);
  cout << " input root = " << inputroot << endl;
  TObjArray HList(0);
  
  TString YtarDeltaCutFile;
  TFile *fYtarDeltaCut;
  vector <TCutG*> ytar_delta_cut;
  if (CutYtarFlag) {
    YtarDeltaCutFile=Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
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
  }
  //
  TString outCutFile;
  TFile *fcut;
  vector<vector<vector<TCutG*> > > ypfp_yfp_cut;
  vector<vector<vector<Int_t> > > ypfp_yfp_cut_flag;
  ypfp_yfp_cut.resize(NumFoil);
  ypfp_yfp_cut_flag.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    ypfp_yfp_cut[nf].resize(ndelcut);
    ypfp_yfp_cut_flag[nf].resize(ndelcut);
  }
  if (CutYpFpYFpFlag) {
    outCutFile=Form("cuts/YpFpYFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    fcut = new TFile(outCutFile);
    cout << " Cut file = " << outCutFile << endl;
    fcut->cd();
    for  (Int_t nf=0;nf<NumFoil;nf++) {
      for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<11;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  } else {
	    //cout << " No hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  }
	}}}
  }
  //
//
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
  if (CutXpFpXFpFlag) {
    xpfp_xfp_outCutFile=Form("cuts/XpFpXFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    xpfp_xfp_fcut = new TFile(xpfp_xfp_outCutFile);
    cout << "xpfp_xfp_ Cut file = " << xpfp_xfp_outCutFile << endl;
    xpfp_xfp_fcut->cd();
    for  (Int_t nf=0;nf<NumFoil;nf++) {
      for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<13;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	    xpfp_xfp_cut[nf][nd].push_back(tempg);
	  } else {
	    //cout << " No hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	    xpfp_xfp_cut[nf][nd].push_back(tempg);
	  }
	}}}
  }
  //
  TFile *fsimc = new TFile(inputroot); 
  TTree *tsimc = (TTree*) fsimc->Get("T");
  // Define branches
  int NMAX = 100000; 

  double yfp[NMAX];
  tsimc->SetBranchAddress("bb.tr.r_y",yfp);
  //tsimc->SetBranchAddress("bb.gem.track.y",yfp);
  double ypfp[NMAX];
  tsimc->SetBranchAddress("bb.tr.r_ph",ypfp);
  //tsimc->SetBranchAddress("bb.gem.track.ph",ypfp);
  double xfp[NMAX];
  tsimc->SetBranchAddress("bb.tr.r_x",xfp);
  //tsimc->SetBranchAddress("bb.gem.track.x",xfp);
  double xpfp[NMAX];
  tsimc->SetBranchAddress("bb.tr.r_th",xpfp);
  //tsimc->SetBranchAddress("bb.gem.track.th",xpfp);
  double ytar_og[NMAX];
  tsimc->SetBranchAddress("bb.tr.tg_y",ytar_og);
  double yptar_og[NMAX];
  tsimc->SetBranchAddress("bb.tr.tg_ph",yptar_og);
  double xtar_og[NMAX];
  tsimc->SetBranchAddress("bb.tg.tg_x",xtar_og);
  double xptar_og[NMAX];
  tsimc->SetBranchAddress("bb.tr.tg_th",xptar_og);
  double ntracks;
  tsimc->SetBranchAddress("bb.gem.track.ntrack",&ntracks);
  double nhits[NMAX];
  tsimc->SetBranchAddress("bb.gem.track.nhits",nhits);
  double chisq[NMAX];
  tsimc->SetBranchAddress("bb.gem.track.chi2ndf",chisq);
  double esumps = 0;
  tsimc->SetBranchAddress("bb.ps.e_c",&esumps);
  double esumsh = 0;
  tsimc->SetBranchAddress("bb.sh.e_c",&esumsh);
  double epratio = 0;
  tsimc->SetBranchAddress("bb.etot_over_p",&epratio);

  //define the variables
  double vx, vy, vz, px, py, pz;
  double p, xptar, yptar, ytar, xtar;
  double p_fit, xptar_fit, yptar_fit, ytar_fit; //Fit is reconstructed using fit coefficients, no smearing for detector resolution
  double pthetabend_fit;
  //double pinv_fit;
  double vz_fit;
  double thetabend_fit;
  double xtar_fit;
  double xsieve, ysieve;
  double z0 = 1.18981;//1.59027;//1.50922;//1.49244;//1.42641;//1.1957;//distance to face of sieve,[m]?, 1.172m is for BB at 1.55m, BB at 1.85m = 1.472m
  TVector3 spec_xaxis_fp,spec_yaxis_fp, spec_zaxis_fp;
  TVector3 spec_xaxis_tgt,spec_yaxis_tgt, spec_zaxis_tgt;
  double tracker_pitch_angle = 0.152622;//0.154066;//0.1514246;//0.15421;//0.15321;//0.152139;//10.22*3.14/180.0;//put this into the input file

  TVector3 BB_zaxis( sin(CentAngle), 0.0, cos(CentAngle) ); //BB is on beam right, global x axis points to beam left
  TVector3 BB_xaxis(0,-1,0); //X axis of transport coordinates is vertically down:
  TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();
  
  spec_xaxis_tgt = BB_xaxis;
  spec_yaxis_tgt = BB_yaxis;
  spec_zaxis_tgt = BB_zaxis;
  
  spec_zaxis_fp = BB_zaxis;
  spec_yaxis_fp = BB_yaxis;
  spec_zaxis_fp.Rotate(-tracker_pitch_angle, spec_yaxis_fp);
  spec_xaxis_fp = spec_yaxis_fp.Cross(spec_zaxis_fp).Unit();
  
  // Define histograms
  cout<<"setup histos"<<endl;
  TH1F *hytar = new TH1F("hytar",Form("Run %d ; Ytar; Counts",nrun),500,-0.25,0.25);
  HList.Add(hytar);
  TH1F *hztar = new TH1F("hztar",Form("Run %d ; Ztar; Counts",nrun),500,-0.3,0.3);
  HList.Add(hztar);
  TH2F *hXptarDelta = new TH2F("hXptarDelta",Form("Run %d ; Xptar ; P x #theta_{bend}",nrun),120,-.6,.6,100,0,1);
  HList.Add(hXptarDelta);
  TH2F *hYptarDelta = new TH2F("hYptarDelta",Form("Run %d ; Yptar ; P x #theta_{bend}",nrun),120,-.4,.4,100,0,1);
  HList.Add(hYptarDelta);
  TH2F *hYtarDelta = new TH2F("hYtarDelta",Form("Run %d ; Ytar ; P x #theta_{bend}",nrun),100,-0.3,0.3,100,0,1);
  HList.Add(hYtarDelta);
  //
  TH2F *hYpFpYFp_all = new TH2F("hYpFpYFp_all",Form("Run %d ; Ypfp ; Yfp",nrun),100,-.3,.3,100,-0.3,0.3);
  HList.Add(hYpFpYFp_all);
  TH2F *hYFpXFp = new TH2F("hYFpXFp",Form("Run %d ; Yfp ; Xfp",nrun),100,-0.3,0.3,100,-0.7,0.7);
  HList.Add(hYFpXFp);
  TH2F *hXpFpXFp = new TH2F("hXpFpXFp",Form("Run %d ; Xpfp ; Xfp",nrun),100,-.7,.7,100,-0.7,0.7);
  HList.Add(hXpFpXFp);
  TH2F *hYtarYptar = new TH2F("hYtarYptar",Form("Run %d ; Yptar ; Ytar",nrun),100,-.3,.3,100,-0.15,0.15);
  HList.Add(hYtarYptar);
  TH2F *hZtarDelta = new TH2F("hZtarDelta",Form("Run %d ; Ztar ; P x #theta_{bend}",nrun),100,-0.3,0.3,100,0,1);
  HList.Add(hZtarDelta);

  TH1F *h_p = new TH1F("h_p",Form("Run %d ; P, recon",nrun),100,0,10);
  HList.Add(h_p);//p
  TH1F *h_pinvtheta = new TH1F("h_pinvtheta",Form("Run %d ; P x #theta_{bend}",nrun),100,0,1);
  HList.Add(h_pinvtheta);//ptheta
  TH2F *hPinvthetaVx = new TH2F("hPinvthetaVx",Form("Run %d ; P x #theta_{bend}; xfp",nrun),100,0,1,100,-0.7,0.7);
  HList.Add(hPinvthetaVx);//ptheta vs x
  TH2F *hPinvthetaVxtar = new TH2F("hPinvthetaVxtar",Form("Run %d ; P x #theta_{bend}; xTar",nrun),100,0,1,100,-0.05,0.05);
  HList.Add(hPinvthetaVxtar);//ptheta vs x
   TH1F *h_ptheta = new TH1F("h_ptheta",Form("Run %d ; P x #theta_{bend}",nrun),100,0,1);
  HList.Add(h_ptheta);//ptheta
  TH2F *hPthetaVx = new TH2F("hPthetaVx",Form("Run %d ; P x #theta_{bend}; xfp",nrun),100,0,1,100,-0.7,0.7);
  HList.Add(hPthetaVx);//ptheta vs x
  TH2F *hPthetaVxtar = new TH2F("hPthetaVxtar",Form("Run %d ; P x #theta_{bend}; xTar",nrun),100,0,1,100,-0.05,0.05);
  HList.Add(hPthetaVxtar);//ptheta vs x
  TH2F *hthetaVp = new TH2F("hthetaVp",Form("Run %d ; #theta_{bend}; P [GeV/c]",nrun),100,0,1,100,0,10);
  HList.Add(hthetaVp);
  TH2F *h_xsVys = new TH2F("h_xsVys",Form("Run %d ; y_{sieve}; x_{sieve}",nrun),100,-0.2,0.2,100,-0.4,0.4);
  HList.Add(h_xsVys);
  TH2F *h_xsVysW = new TH2F("h_xsVysW",Form("Run %d (physics weighted [s^{-1}]) ; y_{sieve}; x_{sieve}",nrun),100,-0.2,0.2,100,-0.4,0.4);
  HList.Add(h_xsVysW);
  
   TH2F *hYpFpYFp_cut0 = new TH2F("hYpFpYFp_cut0",Form("Run %d, yS=2 ; Ypfp ; Yfp",nrun),100,-.3,.3,100,-0.3,0.3);
  HList.Add(hYpFpYFp_cut0);
  TH2F *hXpFpXFp_cut0 = new TH2F("hXpFpXFp_cut0",Form("Run %d, xS=2 ; Xpfp ; Xfp",nrun),100,-.7,.7,100,-0.7,0.7);
  HList.Add(hXpFpXFp_cut0);
  TH2F *h_zp = new TH2F("h_zp",Form("Run %d;zTar;P, recon",nrun),100,-0.3,0.3,100,0,10);
  HList.Add(h_zp);
  TH1F *h_esumNorm = new TH1F("h_esumNorm",Form("Run %d;Esum/P",nrun),100,-2,2);
  HList.Add(h_esumNorm);
  TH1F *h_pthetabend = new TH1F("h_pthetabend",Form("Run %d;p*thetabend",nrun),100,0,1);
  HList.Add(h_pthetabend);
  TH1F *h_xptar = new TH1F("h_xptar",Form("Run %d;xptar",nrun),100,-0.5,0.5);
  HList.Add(h_xptar);
  TH1F *h_yptar = new TH1F("h_yptar",Form("Run %d;yptar",nrun),100,-0.2,0.2);
  HList.Add(h_yptar);

  //
  vector <TH2F*> hYsDelta;
  hYsDelta.resize(NumFoil);
  vector <TH2F*> hXsDelta;
  hXsDelta.resize(NumFoil);
  vector <TH2F*> hYpFpYFp;
  hYpFpYFp.resize(NumFoil);
  vector<vector<vector<TH2F*> > > hYsXs_DelCut_YpYfpCut;
  vector<vector<vector<TH2F*> > > hYsXs_DelCut_XpXfpCut;
  vector<vector<vector<TH1F*> > > hXs_DelCut_YpYfpCut;
  vector<vector<TH2F*> > hYsXs_DelCut;
  vector<vector<TH2F*> > hYpFpYFp_DelCut;
  vector<vector<TH2F*> > hXpFpXFp_DelCut;
  cout << " setup DelCut 2d" << endl;
  hYsXs_DelCut.resize(NumFoil);
  hYsXs_DelCut_YpYfpCut.resize(NumFoil);
  hYsXs_DelCut_XpXfpCut.resize(NumFoil);
  hXs_DelCut_YpYfpCut.resize(NumFoil);
  hYpFpYFp_DelCut.resize(NumFoil);
  hXpFpXFp_DelCut.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    hYsXs_DelCut[nf].resize(ndelcut);
    hYsXs_DelCut_YpYfpCut[nf].resize(ndelcut);
    hYsXs_DelCut_XpXfpCut[nf].resize(ndelcut);
    hXs_DelCut_YpYfpCut[nf].resize(ndelcut);
    hYpFpYFp_DelCut[nf].resize(ndelcut);
    hXpFpXFp_DelCut[nf].resize(ndelcut);
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      hYsXs_DelCut_YpYfpCut[nf][nd].resize(13);
      hYsXs_DelCut_XpXfpCut[nf][nd].resize(13);
      hXs_DelCut_YpYfpCut[nf][nd].resize(13);
    }
  }
  cout << " finish setup Cut 2d" << endl;
  for  (Int_t nc=0;nc<NumFoil;nc++) {
    hYsDelta[nc] = new TH2F(Form("hYsDelta_Foil_%d",nc),Form("Run %d Foil %d; Ys ; P x #theta_{bend}",nc,nrun),100,-0.2,0.2,50,0.,1);
    HList.Add(hYsDelta[nc]);
    hXsDelta[nc] = new TH2F(Form("hXsDelta_Foil_%d",nc),Form("Run %d Foil %d; Xs ; P x #theta_{bend}",nc,nrun),100,-0.4,0.4,50,0,1);
    HList.Add(hXsDelta[nc]);
    hYpFpYFp[nc] = new TH2F(Form("hYpFpYFp_%d",nc),Form("Run %d Foil %d; Ypfp ; Yfp",nrun,nc),100,-.3,.3,100,-.3,.3);
    HList.Add(hYpFpYFp[nc]);
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      hYsXs_DelCut[nc][nd]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d Cut %3.1f; Ys ; Xs",nrun,nc,delcut[nd]),50,-0.2,0.2,100,-0.4,0.4);
      HList.Add(hYsXs_DelCut[nc][nd]);
      for  (Int_t ny=0;ny<13;ny++) {
	hYsXs_DelCut_YpYfpCut[nc][nd][ny]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d_FpCut_%d",nc,nd,ny),Form("Run %d Foil %d Cut %3.1f Ys=%d; Ys ; Xs",nrun,nc,delcut[nd],ny),100,-0.2,0.2,100,-0.4,0.4);
	HList.Add(hYsXs_DelCut_YpYfpCut[nc][nd][ny]);
	hYsXs_DelCut_XpXfpCut[nc][nd][ny]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d_XFpCut_%d",nc,nd,ny),Form("Run %d Foil %d Cut %3.1f Xs=%d; Ys ; Xs",nrun,nc,delcut[nd],ny),100,-0.2,0.2,100,-0.4,0.4);
	HList.Add(hYsXs_DelCut_XpXfpCut[nc][nd][ny]);
	hXs_DelCut_YpYfpCut[nc][nd][ny]  = new TH1F(Form("hXs_Foil_%d_DelCut_%d_FpCut_%d",nc,nd,ny),Form("Run %d Foil %d Cut %3.1f Ys=%d; Xs",nrun,nc,delcut[nd],ny),100,-0.4,0.4);
	HList.Add(hXs_DelCut_YpYfpCut[nc][nd][ny]);
      }
      hYpFpYFp_DelCut[nc][nd]  = new TH2F(Form("hYpFpYFp_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d Cut %3.1f; Ypfp ; Yfp",nrun,nc,delcut[nd]),75,-.3,.3,150,-0.3,0.3);
      HList.Add(hYpFpYFp_DelCut[nc][nd]);
      hXpFpXFp_DelCut[nc][nd]= new TH2F(Form("hXpFpXFp_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d Cut %3.1f; Xpfp ; Xfp",nrun,nc,delcut[nd]),150,-0.7,0.7,150,-0.7,0.7);
      HList.Add(hXpFpXFp_DelCut[nc][nd]);
    }
  }	  
  

  //reading the model file and storing the data in a matrix, M
  //string  modelfilename = "optics_sbs9.txt";
  string  modelfilename = "newfit_sbs9.dat";
  ifstream modelfile(modelfilename.c_str());
  TString currentline;
  //while( currentline.ReadLine(inputfile) ){}
  
  int row_M = 0, col_M = 9;
  modelfile >> row_M;
  TMatrixD M(row_M,col_M);
  for(int row=0; row<row_M; row++){
    for(int col=0; col<col_M; col++){ 
      modelfile >> M(row,col);
      cout<<M(row,col)<<" ";
    }
    cout<<endl;
  }

  // loop over entries
  Long64_t nentries = tsimc->GetEntries();
  cout << " start loop " << nentries << endl;
    for (int i = 0; i < nentries; i++) {
      //  for (int i = 0; i < 1000000; i++) {
    tsimc->GetEntry(i);
    if (i%100000==0)cout << " Entry = " << i << endl;

    //determine if good track
    bool goodtrack = false;
    int itrack = 0;
  
    for (int ii=0; ii<ntracks; ii++){
      //if (nhits[ii]>=4 && xfp[ii]<0.55 && xfp[ii]>-0.55 && chisq[ii]<20.0 && esumps>0.2 && epratio>0.5){
      if (nhits[ii]>=4 && xfp[ii]<0.55 && xfp[ii]>-0.55 && chisq[ii]<30.0 && esumps>0.2 && epratio>0.7&&esumps+esumsh>0.3){
	goodtrack=true;
	itrack = ii;
      }
    }


    vy=0;
    
    if (goodtrack){
      //reconstruct the target quantities
	xtar_fit = -vy;
	
	for( int iter=0; iter<3; iter++ ){

	  xptar_fit = 0.0;
	  yptar_fit = 0.0;
	  ytar_fit = 0.0;
	  pthetabend_fit = 0.0;
	
	  for (int row=0; row<row_M; row++){
	    xptar_fit += M(row,0)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_fit,M(row,8));
	    yptar_fit += M(row,1)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_fit,M(row,8));
	    ytar_fit += M(row,2)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_fit,M(row,8));
	    //pinv_fit += M(row,3)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_fit,M(row,8));
	    pthetabend_fit += M(row,3)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_fit,M(row,8));
	    	
	  }

	  //BB, beam left:
	  vz_fit = -ytar_fit / (sin(CentAngle) + cos(CentAngle)*yptar_fit);	   
	  xtar_fit   = -vy - vz_fit * cos(CentAngle) * xptar_fit;
	}

	//calculate theta bend:
	TVector3 phat_tgt_recon(xptar_fit, yptar_fit, 1.0 );
	phat_tgt_recon = phat_tgt_recon.Unit();

	TVector3 phat_tgt_recon_global = phat_tgt_recon.X() * spec_xaxis_tgt +
	  phat_tgt_recon.Y() * spec_yaxis_tgt +
	  phat_tgt_recon.Z() * spec_zaxis_tgt;

	TVector3 phat_fp_recon(xpfp[itrack], ypfp[itrack], 1.0 );
	phat_fp_recon = phat_fp_recon.Unit();
	
	TVector3 phat_fp_recon_global = phat_fp_recon.X() * spec_xaxis_fp +
	  phat_fp_recon.Y() * spec_yaxis_fp +
	  phat_fp_recon.Z() * spec_zaxis_fp;

	double thetabend_recon = acos( phat_fp_recon_global.Dot( phat_tgt_recon_global ) );

	int pexpansion_flag = 0;
	double p_recon;
	if( pexpansion_flag == 0 ){
	  p_recon = pthetabend_fit/thetabend_recon;
	} else {
	  p_recon = 1.0/pthetabend_fit;
	}
	double pinv_recon = 1.0/p_recon;
	
	xsieve = xtar_fit + xptar_fit*z0;
	ysieve = ytar_fit + yptar_fit*z0;


	//h_esumNorm->Fill((esumps+esumsh)/p_recon);
	h_xptar->Fill(xptar_fit);
	h_yptar->Fill(yptar_fit);

	hytar->Fill(ytar_fit);
	hztar->Fill(vz_fit);
	hXptarDelta->Fill(xptar_fit,pthetabend_fit);
	hYptarDelta->Fill(yptar_fit,pthetabend_fit);
	hYtarDelta->Fill(ytar_fit,pthetabend_fit);
	hYtarYptar->Fill(yptar_fit,ytar_fit);
	hYpFpYFp_all->Fill(ypfp[itrack],yfp[itrack]);
	hXpFpXFp->Fill(xpfp[itrack],xfp[itrack]);
	hYFpXFp->Fill(yfp[itrack],xfp[itrack]); 
	hYtarYptar->Fill(yptar_fit,ytar_fit);
	hZtarDelta->Fill(vz_fit,pthetabend_fit);
	h_zp->Fill(vz_fit,p_recon);
	h_pthetabend->Fill(pthetabend_fit);
	  
	h_p->Fill(p_recon);
	h_pinvtheta->Fill(pthetabend_fit);
	hPinvthetaVx->Fill(pthetabend_fit,xfp[itrack]);
	hPinvthetaVxtar->Fill(pthetabend_fit,xtar_fit);
	
	h_ptheta->Fill(pthetabend_fit);
	hPthetaVx->Fill(pthetabend_fit,xfp[itrack]);
	hPthetaVxtar->Fill(pthetabend_fit,xtar_fit);
	hthetaVp->Fill(thetabend_recon,p_recon);
	h_xsVys->Fill(ysieve,xsieve);

	for  (UInt_t nc=0;nc<ytar_delta_cut.size();nc++) {
	  if (ytar_delta_cut[nc]->IsInside(ytar_fit,pthetabend_fit))	{ 
	    hYsDelta[nc]->Fill(ysieve,pinv_recon*thetabend_fit);
	    hXsDelta[nc]->Fill(xsieve,pinv_recon*thetabend_fit);
	    hYpFpYFp[nc]->Fill(ypfp[itrack],yfp[itrack]);
	    if (nc==0 && abs(ysieve-(2.*0.0381-0.0381*3))<0.03){
	      hYpFpYFp_cut0->Fill(ypfp[itrack],yfp[itrack]);
	    }
	    
	    for  (Int_t nd=0;nd<ndelcut;nd++) {
	      if ( pthetabend_fit >=delcut[nd]-delwidth[nd] && pthetabend_fit <delcut[nd]+delwidth[nd]) {
		hYsXs_DelCut[nc][nd]->Fill(ysieve,xsieve); 
		hYpFpYFp_DelCut[nc][nd]->Fill(ypfp[itrack],yfp[itrack]);
		hXpFpXFp_DelCut[nc][nd]->Fill(xpfp[itrack],xfp[itrack]);
		Int_t f_ny=-1;
		for  (UInt_t ny=0;ny<7;ny++) {
		  if (CutYpFpYFpFlag && ypfp_yfp_cut[nc][nd][ny] && ypfp_yfp_cut[nc][nd][ny]->IsInside(ypfp[itrack],yfp[itrack])) {
		    hYsXs_DelCut_YpYfpCut[nc][nd][ny]->Fill(ysieve,xsieve);
		    hXs_DelCut_YpYfpCut[nc][nd][ny]->Fill(xsieve);
		    f_ny=ny;
		  }
		}
		for  (UInt_t nx=0;nx<13;nx++) {
		  if (f_ny != -1 && CutXpFpXFpFlag && xpfp_xfp_cut[nc][nd][nx] && xpfp_xfp_cut[nc][nd][nx]->IsInside(xpfp[itrack],xfp[itrack])) {
		    hYsXs_DelCut_XpXfpCut[nc][nd][nx]->Fill(ysieve,xsieve);		    
		  }
		}
	      }
	    }
	  }
	}	
    }//end if good track
     
  }//end loop entries

  //
  TFile hsimc(outputhist,"recreate");
  HList.Write();
  //
}

