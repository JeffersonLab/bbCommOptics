#include "TVector3.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TCutG.h"
#include "TEventList.h"
#include "TMinuit.h"
#include <iostream>
#include <fstream>
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
#include <vector>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
using namespace std;

const int nxsieve = 13;
const int nysieve = 7;

const double PI = TMath::Pi();

//Example code to align the GEM stack relative to the spectrometer coordinate system and determine distance from target to sieve plate using
//known dimensions and hole spacings of sieve and known dimensions and internal alignments of GEMs:
//We are going to be lazy and use MINUIT even though this problem could easily be linearized:

const double sigma_pos = 0.0001; //0.1 mm (wild guess)
const double sigma_slope = sigma_pos/1.6; //0.1 mm/1.6 m ~= 6e-5 gives rough order of magnitude for slope uncertainty, of course multiple scattering might invalidate that.

int NTRACKS;

//These arrays will be filled by Holly's code to select tracks going through the individual sieve holes:
vector<double> XTRACK,YTRACK,XPTRACK,YPTRACK; //Measured track parameters in GEM internal coordinate system
vector<double> XSIEVE,YSIEVE; //known x and y center positions of the holes the tracks went through.
//We assume the single foil position defines the global origin and acts as a point source of straight-line rays coming through the sieve holes into the GEMS

void CHI2_FCN( int &npar, double *gin, double &f, double *par, int flag ){
  double chi2 = 0.0;

  //Parameters are:
  
  // zsieve: distance from target to sieve along spectrometer central ray
  // (x0,y0,z0): coordinates of the GEM origin relative to the target center in spectrometer coordinates
  // (theta,phi): polar and azimuthal angles of the Z axis of the GEM stack relative to the spectrometer central ray:
  TVector3 GEMPOS(par[0], par[1], par[2] );
  
  double ZSIEVE = par[3];
  double thetaGEM = par[4];
  double phiGEM = par[5]; 

  //Define unit vectors along the GEM internal axes with respect to the global system:
  TVector3 GEM_zaxis( sin(thetaGEM)*cos(phiGEM), sin(thetaGEM)*sin(phiGEM), cos(thetaGEM) );

  TVector3 Global_yaxis(0,1,0);
  //Not clear if we will need these:
  TVector3 Global_xaxis(1,0,0);
  TVector3 Global_zaxis(0,0,1); 
  
  TVector3 GEM_xaxis = (Global_yaxis.Cross( GEM_zaxis )).Unit();
  TVector3 GEM_yaxis = (GEM_zaxis.Cross( GEM_xaxis )).Unit();
  
  chi2 = 0.0;

  for( int i=0; i<NTRACKS; i++ ){
    //Track position at sieve plate: XSIEVE and YSIEVE have to be determined by Holly's code:
    TVector3 TrackSievePos( XSIEVE[i], YSIEVE[i], ZSIEVE );

    //Track direction 
    TVector3 TrackDirGlobal = TrackSievePos.Unit();

    //Now compute track intersection with plane of first GEM:
    //Since all tracks start at origin, 
    //Equation is: ( s * TrackDir - GEMPOS ) dot GEM_zaxis = 0
    // OR: s = (GEMPOS dot GEM_zaxis)/(TrackDir dot GEM_zaxis)
    
    double sintersect = GEMPOS.Dot( GEM_zaxis )/(TrackDirGlobal.Dot( GEM_zaxis ) );

    TVector3 TrackIntersect_FirstGEMplane = sintersect * TrackDirGlobal;

    double xfp_expect = (TrackIntersect_FirstGEMplane - GEMPOS).Dot( GEM_xaxis );
    double yfp_expect = (TrackIntersect_FirstGEMplane - GEMPOS).Dot( GEM_yaxis );
    double xpfp_expect = TrackDirGlobal.Dot( GEM_xaxis )/TrackDirGlobal.Dot( GEM_zaxis );
    double ypfp_expect = TrackDirGlobal.Dot( GEM_yaxis )/TrackDirGlobal.Dot( GEM_zaxis );

    chi2 += ( pow( xfp_expect - XTRACK[i], 2 ) + pow( yfp_expect - YTRACK[i], 2 ) ) / (sigma_pos*sigma_pos) +
      ( pow( xpfp_expect - XPTRACK[i], 2 ) + pow( ypfp_expect - YPTRACK[i], 2 ) ) / (sigma_slope*sigma_slope);
    
  }
  
  f = chi2;
}

void AlignZeroField( const char *configfilename ){

  TChain *C = new TChain("T");
  
  ifstream configfile(configfilename);

  TString prefix = "bb.gem";

  //These are in spectrometer transport coordinates for an assumed magnet distance of 1.85 m
  double GEMX0=-0.1637, GEMY0=0.0, GEMZ0=2.958;
  double ZSIEVE=1.85-13.37*2.54/100.0; // = 1.51 m
  double GEMtheta = 10.0*PI/180.0; //this will become the pitch angle
  double GEMphi = 180.0*PI/180.0; 

  TCut globalcut = "";
  
  if( configfile ){
    TString currentline;
    
    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	C->Add(currentline.Data());
      }
    }

    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig")){
      if( !currentline.BeginsWith("#") ){
	TObjArray *tokens = currentline.Tokenize(" ");
	
	int ntokens = tokens->GetEntries();

	if( ntokens >= 2 ){
	  TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	  if( skey == "prefix" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    prefix = stemp;
	  }

	  if( skey == "GEMX0" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMX0 = stemp.Atof();
	  }
	  if( skey == "GEMY0" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMY0 = stemp.Atof();
	  }

	  if( skey == "GEMZ0" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMZ0 = stemp.Atof();
	  }

	  if( skey == "ZSIEVE" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    ZSIEVE = stemp.Atof();
	  }

	  if( skey == "GEMtheta" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMtheta = stemp.Atof()*PI/180.0;
	  }

	  if( skey == "GEMphi" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMphi = stemp.Atof()*PI/180.0;
	  }
	  
	}
      }
    }

    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
      if( !currentline.BeginsWith("#") ){
	globalcut += currentline;
      }
    }
  } else {
    return;
  }

  TEventList *elist = new TEventList("elist");
  
  C->Draw(">>elist",globalcut);
  
  cout << "Number of events passing global cut = " << elist->GetN() << endl;

  //TODO: need to declare variables and set tree branch addresses

  UInt_t MAXNTRACKS = 10000;
  
  double ntracks, besttrack;

  //We'll use the GEM track variables directly rather than the BigBite ones to avoid confusion initially.
  //The GEM track variables are guaranteed to be given in the internal GEM coordinates:
  vector<double> tracknhits(MAXNTRACKS);
  vector<double> trackChi2NDF(MAXNTRACKS);
  vector<double> trackX(MAXNTRACKS);
  vector<double> trackY(MAXNTRACKS);
  vector<double> trackXp(MAXNTRACKS);
  vector<double> trackYp(MAXNTRACKS);

  C->SetBranchStatus("*",0);

  TString branchname;
  
  C->SetBranchStatus( branchname.Format("%s.track.ntrack",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.chi2ndf",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.x",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.y",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.xp",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.yp",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.nhits",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.besttrack",prefix.Data() ), 1 );
  //main event loop here:
  C->SetBranchAddress( branchname.Format("%s.track.ntrack",prefix.Data() ), &ntracks );
  C->SetBranchAddress( branchname.Format("%s.track.chi2ndf",prefix.Data() ), &(trackChi2NDF[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.x",prefix.Data() ), &(trackX[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.y",prefix.Data() ), &(trackY[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.xp",prefix.Data() ), &(trackXp[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.yp",prefix.Data() ), &(trackYp[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.nhits",prefix.Data() ), &(tracknhits[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.besttrack",prefix.Data() ), &besttrack );
  
  TObjArray HList(0);
  TH2D *hxyfp = new TH2D("hxyfp"," ; y_{fp} (m) ; x_{fp} (m)",100,-0.3,0.3,250,-1.1,1.1); HList.Add(hxyfp);
  TH2D *hxysieve = new TH2D("hxysieve", " ; y_{sieve} (m) ; x_{sieve} (m)", 200,-0.4,0.4,200,-0.4,0.4); HList.Add(hxysieve);
  TH2F *hYpFpYFp = new TH2F("hYpFpYFp"," ; Ypfp ; Yfp",100,-.3,.3,100,-0.3,0.3); HList.Add(hYpFpYFp);
  TH2F *hXpFpXFp = new TH2F("hXpFpXFp"," ; Xpfp ; Xfp",100,-.7,.7,100,-0.7,0.7); HList.Add(hXpFpXFp);
  TH2D *hxyfp_sel = new TH2D("hxyfp_sel","with sieve cuts; y_{fp} (m) ; x_{fp} (m)",100,-0.3,0.3,250,-1.1,1.1); HList.Add(hxyfp_sel);
  TH2D *hxysieve_sel = new TH2D("hxysieve_sel", "with sieve cuts; y_{sieve} (m) ; x_{sieve} (m)", 200,-0.4,0.4,200,-0.4,0.4); HList.Add(hxysieve_sel);
  TH2F *hYpFpYFp_sel = new TH2F("hYpFpYFp_sel","with sieve cuts; Ypfp ; Yfp",100,-.3,.3,100,-0.3,0.3); HList.Add(hYpFpYFp_sel);
  TH2F *hXpFpXFp_sel = new TH2F("hXpFpXFp_sel","with sieve cuts; Xpfp ; Xfp",100,-.7,.7,100,-0.7,0.7); HList.Add(hXpFpXFp_sel);
  TH2D *hxysieve_new = new TH2D("hxysieve_new", "post fit; y_{sieve} (m) ; x_{sieve} (m)", 200,-0.4,0.4,200,-0.4,0.4); HList.Add(hxysieve_new);

  long nevent=0;

  TVector3 Global_yaxis(0,1,0);
  //Not clear if we will need these:
  TVector3 Global_xaxis(1,0,0);
  TVector3 Global_zaxis(0,0,1);

  TVector3 GEMzaxis( sin(GEMtheta)*cos(GEMphi), sin(GEMtheta)*sin(GEMphi), cos(GEMtheta) );
  TVector3 GEMxaxis = (Global_yaxis.Cross(GEMzaxis)).Unit();
  TVector3 GEMyaxis = (GEMzaxis.Cross(GEMxaxis)).Unit();
  TVector3 GEMorigin( GEMX0, GEMY0, GEMZ0 );
  TVector3 SieveOrigin(0,0,ZSIEVE);

  TString outCutFile;
  outCutFile=Form("sieve_cut.root");
  cout << "outCutFile =  " << outCutFile << endl;
  TFile fcut(outCutFile,"UPDATE");
  TString hname_cut;

  //Here are the desired sieve positions:
  vector <Double_t> xs_cent{-(0.3+0.0492)+0.0493/cos(18.*PI/180.),
      -(0.3+0.0492)+(0.0493+0.0492)/cos(18.*PI/180),
      -(0.3+0.0492)+0.1493/cos(9.*PI/180.),
      -(0.3+0.0492)+(0.1493+0.0492)/cos(9.*PI/180.),
      -(0.3+0.0492)+(0.1493+0.0492*2.)/cos(9.*PI/180.),
      -0.0492,
      0.0,
      0.0492,
      0.3+0.0492-(0.1493+0.0492*2.)/cos(9.*PI/180.),
      0.3+0.0492-(0.1493+0.0492)/cos(9.*PI/180.),
      0.3+0.0492-0.1493/cos(9.*PI/180.),
      0.3+0.0492-(0.0493+0.0492)/cos(18.*PI/180),
      0.3+0.0492-0.0493/cos(18.*PI/180.)};
   
  vector <Double_t> ys_cent;
  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=nys*0.0381-0.0381*3;//old sieve
    ys_cent.push_back(pos);
  }
  
  while( C->GetEntry( elist->GetEntry( nevent++ ) ) ){
    
    //do stuff: grab track info from tree, apply cuts, fill the track arrays defined near the top of this macro (see XTRACK, YTRACK, XPTRACK, YPTRACK, XSIEVE, YSIEVE above)
    int itr = int(besttrack);
    if( besttrack >= 0 && ntracks <= MAXNTRACKS ){
      hxyfp->Fill( trackY[itr], trackX[itr] );
      hYpFpYFp->Fill(trackYp[itr],trackY[itr]);
      hXpFpXFp->Fill(trackXp[itr],trackX[itr]);
  
      TVector3 TrackPos_local( trackX[itr], trackY[itr], 0.0 );
      TVector3 TrackDir_local( trackXp[itr], trackYp[itr], 1.0 );
      TrackDir_local = TrackDir_local.Unit();
      
      TVector3 TrackDir_global = TrackDir_local.X() * GEMxaxis +
	TrackDir_local.Y() * GEMyaxis +
	TrackDir_local.Z() * GEMzaxis;

      TVector3 TrackPos_global = GEMorigin + TrackPos_local.X() * GEMxaxis + TrackPos_local.Y() * GEMyaxis + TrackPos_local.Z() * GEMzaxis; 
      
      //Now compute intersection of track ray in global coordinates with sieve slit:
      // (trackpos + s * trackdir - sieveorigin) dot globalzaxis = 0
      // --> s * trackdir dot globalzaxis = (sieveorigin - trackpos) dot globalzaxis
      double sintersect = (SieveOrigin - TrackPos_global).Dot( Global_zaxis ) /
	( TrackDir_global.Dot( Global_zaxis ) );
      TVector3 TrackSievePos = TrackPos_global + sintersect * TrackDir_global;

      hxysieve->Fill( TrackSievePos.Y(), TrackSievePos.X() );
      //cout<<"ysieve: "<<TrackSievePos.Y()<<" xsieve: "<<TrackSievePos.X()<<endl;
    }
  }//end while event

  //Set the sieve hole selection
  //gROOT->Reset();
  TCanvas *histView_Cut; 
  histView_Cut= new TCanvas("histView_Cut","cut",700,500);
  histView_Cut->Divide(1,1);
  histView_Cut->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogz(); 
  
  TCutG*t;
  histView_Cut->Clear();
  hxysieve->Draw("colz");
  hxysieve->SetMinimum(1.);	
  //hxysieve->GetYaxis()->SetRangeUser(-1,1);
  //hxysieve->GetXaxis()->SetRangeUser(-1,1);
  //draw the ideal positions
  TLine* ys_line[nysieve];
  TText* ys_text[nysieve];
  TLine* xs_line[nxsieve];
  TText* xs_text[nxsieve];
  
  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=ys_cent[nys];//for old nysieve=7	
    ys_line[nys]= new TLine(pos,-0.4,pos,0.4);
    ys_text[nys]= new TText(pos,-0.5,Form("%d",nys));
    ys_text[nys]->SetTextColor(2);
    ys_line[nys]->SetLineColor(2);
    ys_line[nys]->SetLineWidth(1);
    ys_line[nys]->Draw("same");
    ys_text[nys]->Draw("same");
  }
  for (Int_t nxs=0;nxs<nxsieve;nxs++) {
    Double_t pos=xs_cent[nxs];//for old nysieve=7	
    xs_line[nxs]= new TLine(-0.4,pos,0.4,pos);
    xs_text[nxs]= new TText(-0.5,pos,Form("%d",nxs));
    xs_text[nxs]->SetTextColor(2);
    xs_line[nxs]->SetLineColor(2);
    xs_line[nxs]->SetLineWidth(1);
    xs_line[nxs]->Draw("same");
    xs_text[nxs]->Draw("same");
  }
  
  histView_Cut->Update();

    
  //in some kind of loop
  int nloop=0;
  while (nloop !=-1 ) {
    fcut.cd();
    for  (Int_t nys=0;nys<nysieve;nys++) {
      for (Int_t nxs=0; nxs<nxsieve;nxs++){
	hname_cut= Form("sievecut_%d_%d",nys,nxs);
	t=(TCutG*)gROOT->FindObject(hname_cut);
	if(t) {
	  // t->Draw("same");
	  t->SetLineColor(1);
	  Double_t xcut,ycut;
	  t->GetPoint(0,xcut,ycut);
	  TText* ystext = new TText(xcut,ycut,Form("(%d,%d)",nys,nxs));
	  //ystext->Draw();
	  histView_Cut->Update();
	}
      }
      }
    
    Int_t yscol=0;
    Int_t xscol=0;
    
    cout <<" Action ( 0 ( set cut) , -10 delete cut,  -100 exit) "  << endl;
    cin >> nloop ;
    if (nloop == -100) break;
    if (!(nloop == -10 || nloop==0 || nloop==-1 || nloop==-100)) continue;//return;
    if (nloop==-10 || nloop ==0) {
      cout << " Which xsieve hole ? " << " xscol = " << xscol << endl;
      cin >> xscol;
      cout << " Which ysieve hole ? " << " yscol = " << yscol << endl;
      cin >> yscol;
      if (yscol >=nysieve)yscol=0;
      if (xscol >=nxsieve)xscol=0;
      hname_cut= Form("sievecut_%d_%d",yscol,xscol);
    }
    if(nloop==-10) {
      fcut.cd();
      fcut.Delete(hname_cut+";*");
      fcut.Delete(hname_cut);
      t=(TCutG*)gROOT->FindObject(hname_cut);
      gROOT->Remove(t);
      t=(TCutG*)gROOT->FindObject(hname_cut);
      if (!t) cout << " delete cut = " <<hname_cut  << endl;
      if (t) cout << " delete cut? = " <<hname_cut  << endl;
      fcut.Write("",TObject::kOverwrite);
    }
    if (nloop==0) {
      t=(TCutG*)gROOT->FindObject(hname_cut);	      
      if (t) {
	fcut.cd();
	fcut.Delete(hname_cut);
	hname_cut= Form("sievecut_%d_%d",yscol,xscol);
	t=(TCutG*)gROOT->FindObject(hname_cut);
	gROOT->Remove(t);
      }
      TCutG*cutg=(TCutG*)gPad->WaitPrimitive("CUTG","CutG");
      histView_Cut->Update();
      TCutG*tmpg= (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
      hname_cut= Form("sievecut_%d_%d",yscol,xscol);
      TCutG*mycutg=(TCutG*)(tmpg->Clone(hname_cut));
      fcut.cd();
      mycutg->Write("",TObject::kOverwrite);
 
    }
    
  }//end while for setting cuts
  //TFile hsimc(outCutFile,"recreate");
  //HList.Write();
  fcut.Close();
    
  //load up the sieve hole cuts
  nevent = 0;
  vector<vector<TCutG*>> sieve_cut;
  sieve_cut.resize(nysieve);
  
  TFile *fcutfileroot;
  fcutfileroot = new TFile("sieve_cut.root");
  fcutfileroot->cd();
  for (Int_t iy=0; iy<nysieve; iy++){
    for (Int_t ix=0;ix<nxsieve;ix++){
      TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("sievecut_%d_%d;1",iy,ix));
      if (tempg)  {
	Int_t npt = tempg->GetN();
	sieve_cut[iy].push_back(tempg);
      }
      else{
	sieve_cut[iy].push_back(tempg);
      }
    }} 
   
  
  while( C->GetEntry( elist->GetEntry( nevent++ ) ) ){
    //In this second loop, make explicit the association between sieve holes and tracks:

    int itr = int(besttrack);
    if( besttrack >= 0 && ntracks <= MAXNTRACKS ){
      
      TVector3 TrackPos_local( trackX[itr], trackY[itr], 0.0 );
      TVector3 TrackDir_local( trackXp[itr], trackYp[itr], 1.0 );
      TrackDir_local = TrackDir_local.Unit();
      
      TVector3 TrackDir_global = TrackDir_local.X() * GEMxaxis +
	TrackDir_local.Y() * GEMyaxis +
	TrackDir_local.Z() * GEMzaxis;

      TVector3 TrackPos_global = GEMorigin + TrackPos_local.X() * GEMxaxis + TrackPos_local.Y() * GEMyaxis + TrackPos_local.Z() * GEMzaxis; 
      
      //Now compute intersection of track ray in global coordinates with sieve slit:
      // (trackpos + s * trackdir - sieveorigin) dot globalzaxis = 0
      // --> s * trackdir dot globalzaxis = (sieveorigin - trackpos) dot globalzaxis
      double sintersect = (SieveOrigin - TrackPos_global).Dot( Global_zaxis ) /
	( TrackDir_global.Dot( Global_zaxis ) );
      TVector3 TrackSievePos = TrackPos_global + sintersect * TrackDir_global;

      int iy_found=-1;
      int ix_found=-1;
      
      for (Int_t iy=0; iy<nysieve; iy++){
	for (Int_t ix=0;ix<nxsieve;ix++){
	  if (sieve_cut[iy][ix]){
	    if (sieve_cut[iy][ix]->IsInside(TrackSievePos.Y(),TrackSievePos.X())){
	      iy_found = iy; ix_found=ix;
	      hxysieve_sel->Fill( TrackSievePos.Y(), TrackSievePos.X() );
	    }
	  }
	}
      }
      //we have an event going through a hole
      if (iy_found != -1 && ix_found != -1){
	hxyfp_sel->Fill( trackY[itr], trackX[itr] );
	hYpFpYFp_sel->Fill(trackYp[itr],trackY[itr]);
	hXpFpXFp_sel->Fill(trackXp[itr],trackX[itr]);
	XTRACK.push_back(trackX[itr]);
	YTRACK.push_back(trackY[itr]);
	XPTRACK.push_back(trackXp[itr]);
	YPTRACK.push_back(trackYp[itr]);
	XSIEVE.push_back(xs_cent[ix_found]);
	YSIEVE.push_back(ys_cent[iy_found]);	
      }
    }
  }
  NTRACKS=XTRACK.size();
  cout<<"found "<<NTRACKS<<" events."<<endl;
  
  //do the thing
  TMinuit *FitZeroField = new TMinuit( 6 );
  
  FitZeroField->SetFCN( CHI2_FCN );
  
  int ierflg = 0;
  
  FitZeroField->mnparm( 0, "GEMX0", GEMX0, 0.1,0,0,ierflg ); //guesstimate 1 cm as initial position accuracy
  FitZeroField->mnparm( 1, "GEMY0", GEMY0, 0.03,0,0,ierflg ); //start with 3-cm uncertainty for y position
  FitZeroField->mnparm( 2, "GEMZ0", GEMZ0, 0.01,0,0,ierflg );
  FitZeroField->mnparm( 3, "ZSIEVE", ZSIEVE, 0.01,0,0,ierflg );
  FitZeroField->mnparm( 4, "GEMtheta", GEMtheta, 1.*PI/180.0, 0, 0, ierflg ); //guesstimate 0.3 degrees as initial angular accuracy
  FitZeroField->mnparm( 5, "GEMphi", GEMphi, 0.3*PI/180.0, 0, 0, ierflg );
    
  double arglist[10];
  arglist[0]=1;
  FitZeroField->mnexcm("SET ERR",arglist,1,ierflg);
  
  arglist[0] = 5000;
  arglist[1] = 1.;

  cout<<"filled all variables for the minimization!"<<endl;
  
  //Uncomment this line when we are actually set up to do a fit:   
  FitZeroField->mnexcm("MIGRAD",arglist,2,ierflg);
  
  //TODO: grab parameters, write them out to file. Profit.
  double GEMX0_fit=-0.205104;//-0.213869;//-0.206968;//-0.198518;//-0.211922;//-0.212184;
  double GEMY0_fit=0.0423978;//0.0394128;//0.0412156;//0.0425868;//0.0406949;//0.0322128;
  double GEMZ0_fit=2.70631;//3.12109;//3.0399;//2.71764;//3.01701;//3.00532;
  double ZSIEVE_fit=1.18981;//1.59027;//1.50922;//1.22518;//1.49244;//1.47786;
  double GEMtheta_fit=0.152622;//0.154066;//0.151425;//0.149694;//0.154208;//0.153710;
  double GEMphi_fit=3.12074;//3.13248;//3.15825;//3.11352;//3.13954;//3.14944;

  //check the results:
  TVector3 GEMzaxisFit( sin(GEMtheta_fit)*cos(GEMphi_fit), sin(GEMtheta_fit)*sin(GEMphi_fit), cos(GEMtheta_fit) );
  TVector3 GEMxaxisFit = (Global_yaxis.Cross(GEMzaxis)).Unit();
  TVector3 GEMyaxisFit = (GEMzaxis.Cross(GEMxaxis)).Unit();
  TVector3 GEMoriginFit( GEMX0_fit, GEMY0_fit, GEMZ0_fit );
  TVector3 SieveOriginFit(0,0,ZSIEVE_fit);

  cout<<"got new stuff"<<endl;
  nevent=0;
  while( C->GetEntry( elist->GetEntry( nevent++ ) ) ){
    //do stuff: grab track info from tree, apply cuts, fill the track arrays defined near the top of this macro (see XTRACK, YTRACK, XPTRACK, YPTRACK, XSIEVE, YSIEVE above)
    int itr = int(besttrack);
    if( besttrack >= 0 && ntracks <= MAXNTRACKS ){
  
      TVector3 TrackPos_local( trackX[itr], trackY[itr], 0.0 );
      TVector3 TrackDir_local( trackXp[itr], trackYp[itr], 1.0 );
      TrackDir_local = TrackDir_local.Unit();
      
      TVector3 TrackDir_global = TrackDir_local.X() * GEMxaxisFit +
	TrackDir_local.Y() * GEMyaxisFit +
	TrackDir_local.Z() * GEMzaxisFit;

      TVector3 TrackPos_global = GEMoriginFit + TrackPos_local.X() * GEMxaxisFit + TrackPos_local.Y() * GEMyaxisFit + TrackPos_local.Z() * GEMzaxisFit; 
      
      //Now compute intersection of track ray in global coordinates with sieve slit:
      // (trackpos + s * trackdir - sieveorigin) dot globalzaxis = 0
      // --> s * trackdir dot globalzaxis = (sieveorigin - trackpos) dot globalzaxis
      double sintersect = (SieveOriginFit - TrackPos_global).Dot( Global_zaxis ) /
	( TrackDir_global.Dot( Global_zaxis ) );
      TVector3 TrackSievePos = TrackPos_global + sintersect * TrackDir_global;

      hxysieve_new->Fill( TrackSievePos.Y(), TrackSievePos.X() );
      //cout<<"ysieve: "<<TrackSievePos.Y()<<" xsieve: "<<TrackSievePos.X()<<endl;
    }
  }//end while event


  TCanvas *histView_fin; 
  histView_fin= new TCanvas("histView_fin","cut",700,500);
  histView_fin->Divide(1,1);
  histView_fin->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogz(); 
  
  histView_fin->Clear();
  hxysieve_new->Draw("colz");
  hxysieve_new->SetMinimum(1.);	

  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=ys_cent[nys];//for old nysieve=7	
    ys_line[nys]= new TLine(pos,-0.4,pos,0.4);
    ys_text[nys]= new TText(pos,-0.5,Form("%d",nys));
    ys_text[nys]->SetTextColor(2);
    ys_line[nys]->SetLineColor(2);
    ys_line[nys]->SetLineWidth(1);
    ys_line[nys]->Draw("same");
    ys_text[nys]->Draw("same");
  }
  for (Int_t nxs=0;nxs<nxsieve;nxs++) {
    Double_t pos=xs_cent[nxs];//for old nysieve=7	
    xs_line[nxs]= new TLine(-0.4,pos,0.4,pos);
    xs_text[nxs]= new TText(-0.5,pos,Form("%d",nxs));
    xs_text[nxs]->SetTextColor(2);
    xs_line[nxs]->SetLineColor(2);
    xs_line[nxs]->SetLineWidth(1);
    xs_line[nxs]->Draw("same");
    xs_text[nxs]->Draw("same");
  }


  
  TFile hsimc("strt_trk_plots.root","recreate");
  HList.Write();




  
 
}//end program
