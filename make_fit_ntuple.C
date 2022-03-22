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

void make_fit_ntuple(Int_t nrun=1814,Int_t FileID=-2){
  Double_t yMP = 0.0;
  Double_t xMP = 0.0;
  const int nxsieve = 13;
  const int nysieve = 7;
  
  Bool_t CutYtarFlag=kTRUE;
  Bool_t CutYpFpYFpFlag=kTRUE;
  Bool_t CutXpFpXFpFlag=kTRUE;
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
  for (Int_t nf=0;nf<NumFoil;nf++) {
    cout << nf << " foil = " << ztar_foil[nf] << endl;
  }
   
  vector <Double_t> ys_cent;
  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=nys*0.0381-0.0381*3;//old sieve
    cout<<"nys: "<<nys<<" pos: "<<pos<<endl;
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

  //
  //
  //
  TString inputroot;
  TString outputroot;
  if (nrun==13674){  
inputroot=Form("Rootfiles/combined_13674_5foil_new.root");//shms_replay_matrixopt_%s_%d.root",OpticsID.Data(),FileID);
  }
  else{inputroot=Form("files/replay_sbs8_13436.root");}
  outputroot= Form("hist/Optics_%s_%d_fit_tree.root",OpticsID.Data(),FileID);
  //
  //
  TString YtarDeltaCutFile;
  TFile *fYtarDeltaCut;
  vector <TCutG*> ytar_delta_cut;
  if (CutYtarFlag) {
    YtarDeltaCutFile=Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
    fYtarDeltaCut = new TFile(YtarDeltaCutFile);
    cout << " Cut file = " << YtarDeltaCutFile << endl;
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
        for (Int_t nc=0;nc<nysieve;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    Int_t npt = tempg->GetN();
	    //cout << "hYpFpYFp_cut = " << nf << " " << nd << " " << nc << " npts = " << npt << endl;
	    //cout << "hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  } else {
	    //cout << " No hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  }
	}}}
  }
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
        for (Int_t nc=0;nc<nxsieve;nc++) {

	  if (nf<1){
	    TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,1,nd));
	    if (tempg)  {
	      cout << "hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	      xpfp_xfp_cut[nf][nd].push_back(tempg);
	    } else {
	      cout << " No hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	      xpfp_xfp_cut[nf][nd].push_back(tempg);
	    }
	  }
	  else{
	    TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	    if (tempg)  {
	      cout << "hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	      xpfp_xfp_cut[nf][nd].push_back(tempg);
	    } else {
	      cout << " No hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	      xpfp_xfp_cut[nf][nd].push_back(tempg);
	    }
	  }
	}//end loop nxsieve
	cout<<"done xsieve loop"<<endl;
      }//end loop ndelcut
      cout<<"done delcut loop"<<endl;
    }//end loop foil
    cout<<"done foils"<<endl;
  }
  cout<<"on to the next great thing"<<endl;
  //
  TFile *fsimc = new TFile(inputroot); 
  TTree *tsimc = (TTree*) fsimc->Get("T");
  // Define branches
  int NMAX = 100000; 

  double yfp[NMAX];
  tsimc->SetBranchAddress("bb.tr.r_y",yfp);
  double ypfp[NMAX];
  tsimc->SetBranchAddress("bb.tr.r_ph",ypfp);
  double xfp[NMAX];
  tsimc->SetBranchAddress("bb.tr.r_x",xfp);
  double xpfp[NMAX];
  tsimc->SetBranchAddress("bb.tr.r_th",xpfp);
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
  tsimc->SetBranchAddress("bb.etot_over_p", &epratio);
  
  //define the variables
  double vx, vy, vz, px, py, pz;
  double p, xptar, yptar, ytar, xtar;
  double p_fit, xptar_fit, yptar_fit, ytar_fit; //Fit is reconstructed using fit coefficients, no smearing for detector resolution
  double p_recon, xptar_recon, yptar_recon, ytar_recon; //recon is reconstructed using fit coefficients, fp quantities smeared by det. resolution
  double pthetabend_true;
  double pthetabend_fit, pthetabend_recon;
  double pinv_fit, pinv_recon;
  double vz_fit, vz_recon;
  double thetabend_true;
  double thetabend_fit;
  double thetabend_recon;
  double xtar_recon, xtar_fit;
  double xsieve, ysieve;
  double z0 = 1.18981;//1.59027;//1.50922;//1.49244;//1.42641;//1.1957;//distance to face of sieve,[m]?
  TVector3 spec_xaxis_fp,spec_yaxis_fp, spec_zaxis_fp;
  TVector3 spec_xaxis_tgt,spec_yaxis_tgt, spec_zaxis_tgt;
  double tracker_pitch_angle = 0.152622;//0.154066;//0.1514246;//0.15421;//0.15321;//0.152139;//10.0*3.14/180.0;//put this into the input file

  Double_t xptarT,ytarT,yptarT,ysieveT,xsieveT,ztarT,ztar,pinvtheta,pmom,weight;
  //
  TFile hroot(outputroot,"recreate");
  TTree *otree = new TTree("TFit","FitTree");
  otree->Branch("ys",&ysieve);
  otree->Branch("ysT",&ysieveT);
  otree->Branch("xs",&xsieve);
  otree->Branch("xsT",&xsieveT);
  otree->Branch("ztarT",&ztarT);
  otree->Branch("xtar",&xtar_fit);
  otree->Branch("ztar",&ztar);
  otree->Branch("xptar",&xptar_fit);
  otree->Branch("yptar",&yptar_fit);
  otree->Branch("ytar",&ytar_fit);
  otree->Branch("xptarT",&xptarT);
  otree->Branch("yptarT",&yptarT);
  otree->Branch("ytarT",&ytarT);
  otree->Branch("xpfp",xpfp);
  otree->Branch("ypfp",ypfp);
  otree->Branch("xfp",xfp);
  otree->Branch("yfp",yfp);
  otree->Branch("pinvtheta",&pinvtheta);
  otree->Branch("pmom",&pmom);
  otree->Branch("weight",&weight);
  otree->Branch("pthetabend",&pthetabend_fit);

  CentAngle=CentAngle*3.14159/180.;
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


  //reading the model file and storing the data in a matrix, M
  string modelfilename = "optics_sbs9.txt";
  //string  modelfilename = "newfit.dat";
  ifstream modelfile(modelfilename.c_str());
  TString currentline;
  //while( currentline.ReadLine(inputfile) ){}
  
  int row_M = 0, col_M = 9;
  modelfile >> row_M;
  TMatrixD M(row_M,col_M);
  
  for(int row=0; row<row_M; row++){
    for(int col=0; col<col_M; col++){ 
      modelfile >> M(row,col);
    }
  }

  // loop over entries
  Double_t zdis_sieve = z0;//front of sieve
  Long64_t nentries = tsimc->GetEntries();
  cout << " start loop " << nentries << endl;
  for (int i = 0; i < nentries; i++) {
    tsimc->GetEntry(i);
    if (i%50000==0)cout << " Entry = " << i << endl;

    //determine if good track
    bool goodtrack = false;
    int itrack = 0;
  
    for (int ii=0; ii<ntracks; ii++){
      if (nhits[ii]>=4 && xfp[ii]<0.55 && xfp[ii]>-0.55 && chisq[ii]<30.0 && esumps>0.2 &&epratio>0.7 && esumps+esumsh>0.3){
	goodtrack=true;
	itrack = ii;
      }
    }
 
    vy = 0;
    if (goodtrack){
      //reconstruct the target quantities
	xtar_fit = -vy;
	xtar_recon = -vy;
	
	for( int iter=0; iter<3; iter++ ){

	  xptar_fit = 0.0;
	  yptar_fit = 0.0;
	  ytar_fit = 0.0;
	  pthetabend_fit = 0.0;
	  pinv_fit = 0.0;
	
	
	  for (int row=0; row<row_M; row++){
	    xptar_fit += M(row,0)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_recon,M(row,8));
	    yptar_fit += M(row,1)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_recon,M(row,8));
	    ytar_fit += M(row,2)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_recon,M(row,8));
	    pinv_fit += M(row,3)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_recon,M(row,8));
	    pthetabend_fit += M(row,3)*pow(xfp[itrack],M(row,4))*pow(yfp[itrack],M(row,5))*pow(xpfp[itrack],M(row,6))*pow(ypfp[itrack],M(row,7))*pow(xtar_recon,M(row,8));
	    
	
	  }

	  //beam left:
	  vz_fit = -ytar_fit / (sin(CentAngle) + cos(CentAngle)*yptar_fit);
	  //beam right:
	  //vz_fit = ytar_fit / (sin(CentAngle) - cos(CentAngle)*yptar_fit);
	  //vz_recon = ytar_recon / (sin(CentAngle) - cos(CentAngle)*yptar_recon);

	  
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

	thetabend_recon = acos( phat_fp_recon_global.Dot( phat_tgt_recon_global ) );

	int pexpansion_flag = 0;
	if( pexpansion_flag == 0 ){
	  p_recon = pthetabend_fit/thetabend_fit;
	} else {
	  p_recon = 1.0/pthetabend_fit;
	}
	pinv_recon = 1.0/p_recon;
	
	xsieve = xtar_fit + xptar_fit*z0;
	ysieve = ytar_fit + yptar_fit*z0;

	Int_t nf_found=-1, nd_found=-1,ny_found=-1,nx_found=-1;
	for  (UInt_t nf=0;nf<ytar_delta_cut.size();nf++) {
	  if (ytar_delta_cut[nf]->IsInside(ytar_fit,pthetabend_fit)) nf_found=nf;
	} 
	for  (UInt_t nd=0;nd<ndelcut;nd++) {
	  if ( pthetabend_fit >=delcut[nd]-delwidth[nd] && pthetabend_fit <delcut[nd]+delwidth[nd])  nd_found=nd;
	}
	if (nf_found!=-1 && nd_found!=-1) {
	  for  (UInt_t ny=0;ny<nysieve;ny++) {
	    if (ypfp_yfp_cut[nf_found][nd_found][ny] && ypfp_yfp_cut[nf_found][nd_found][ny]->IsInside(ypfp[itrack],yfp[itrack])) ny_found=ny;
	  }
	  for  (UInt_t nx=0;nx<nxsieve;nx++) {
	    if (xpfp_xfp_cut[nf_found][nd_found][nx] && xpfp_xfp_cut[nf_found][nd_found][nx]->IsInside(xpfp[itrack],xfp[itrack])) nx_found=nx;
	  }
	}
	if (nf_found !=-1 && nd_found!=-1 && ny_found!=-1 && nx_found!=-1) {
	  //BB beam left:
	  yptarT = (ys_cent[ny_found]+ztar_foil[nf_found]*TMath::Sin(CentAngle))/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
	  xptarT = (xs_cent[nx_found])/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
	  ytarT = -ztar_foil[nf_found]*(TMath::Sin(CentAngle)+yptarT*TMath::Cos(CentAngle));

	  //BB beam left with tgt offsets and mis-pointing:
	  //yptarT = (ys_cent[ny_found]+ztar_foil[nf_found]*TMath::Sin(CentAngle)-reactx*TMath::Cos(CentAngle)+yMP)/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle)-reactx*TMath::Sin(CentAngle));
	  //xptarT = (xs_cent[nx_found]+reacty+xMP)/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle)-reactx*TMath::Sin(CentAngle));
	  //ytarT = -ztar_foil[nf_found]*(TMath::Sin(CentAngle)+yptarT*TMath::Cos(CentAngle))-yptarT*reactx*TMath::Sin(CentAngle))+reactx*TMath::Cos(CentAngle)-yMP;
	  //xtarT = -reacty - xMP - xptarT*ztar_foil[nf_found]*TMath::Cos(CentAngle);`
	  
	  
	  //BB beam right:
	  //yptarT = (ys_cent[ny_found]-ztar_foil[nf_found]*TMath::Sin(CentAngle))/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
	  //xptarT = (xs_cent[nx_found])/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
	  //ytarT = -ztar_foil[nf_found]*(-TMath::Sin(CentAngle)+yptarT*TMath::Cos(CentAngle));


	  pmom=p_fit;
	  ysieveT=ys_cent[ny_found];
	  xsieveT=xs_cent[nx_found];
	  ztarT=ztar_foil[nf_found];
	  ztar=vz_fit;//(ytar-yMP-reactx*(cos(CentAngle)-yptar*sin(CentAngle)))/(-sin(CentAngle)-cos(CentAngle)*yptar)
	  pinvtheta =pthetabend_fit;

	  
	  otree->Fill();
	}
    }//end if good track
  }//end entries
  otree->Write();
}
