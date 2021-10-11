#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2F.h>
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
#include <vector>
using namespace std;


void set_xpfp_xfp_cuts(Int_t nrun=1814,Int_t FileID=-2,Double_t hist_minZ=1.) {
  
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
  const int nxsieve = 13;
  const int nysieve = 7;
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
  CentAngle=CentAngle*3.14159/180.;
  //
  TString inputroot;
  TFile *fhistroot;
  inputroot=Form("hist/Optics_%s_%d_hist.root",OpticsID.Data(),FileID);
  cout << " infile root = " << inputroot << endl;
  fhistroot =  new TFile(inputroot);
  //
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
  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=nys*0.0381-0.0381*3;//old sieve
    ys_cent.push_back(pos);
  }
  
  //
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
  /*
  string oldcoeffsfilename="shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_forward_global.dat";
  ifstream oldcoeffsfile(oldcoeffsfilename.c_str());
  vector<Double_t> xfpcoeffs;
  vector<Double_t> xpfpcoeffs;
  vector<Double_t> yfpcoeffs;
  vector<Double_t> ypfpcoeffs;
  vector<Double_t> lencoeffs;
  vector<Int_t> xtarexpon;
  vector<Int_t> xptarexpon;
  vector<Int_t> ytarexpon;
  vector<Int_t> yptarexpon;
  vector<Int_t> deltaexpon;
  TString currentline;
  int num_recon_terms=0;
  while( currentline.ReadLine(oldcoeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){

    TString sc1(currentline(1,14));
    TString sc2(currentline(15,14));
    TString sc3(currentline(29,14));
    TString sc4(currentline(43,14));
    TString sc5(currentline(57,14));
    
    xfpcoeffs.push_back(sc1.Atof());
    xpfpcoeffs.push_back(sc2.Atof());
    yfpcoeffs.push_back(sc3.Atof());
    ypfpcoeffs.push_back(sc4.Atof());
    lencoeffs.push_back(sc5.Atof());
    
    int expontemp[6];
    
    for(int expon=0; expon<6; expon++){
      TString stemp(currentline(72+expon,1));
      expontemp[expon] = stemp.Atoi();
    }
    
  

    xtarexpon.push_back(expontemp[0]);
    xptarexpon.push_back(expontemp[1]);
    ytarexpon.push_back(expontemp[2]);
    yptarexpon.push_back(expontemp[3]);
    deltaexpon.push_back(expontemp[5]);
    cout << xfpcoeffs[num_recon_terms] << " " <<  xpfpcoeffs[num_recon_terms]<< " " <<  yfpcoeffs[num_recon_terms] << " " <<  ypfpcoeffs[num_recon_terms]<< " " << xtarexpon[num_recon_terms] << xptarexpon[num_recon_terms] << ytarexpon[num_recon_terms] << yptarexpon[num_recon_terms]<< expontemp[4]<< deltaexpon[num_recon_terms]  << endl;
   
    
    num_recon_terms++;
  }
  //
  */
 
  vector<vector<TH2F*>> hXpFpXFp_DelCut;
  vector<TH2F*> temp2d;
  
  TH2F* th;
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
  
  TString hname_cut;
  hname_cut= Form("hXpFpXFp_cut");
    

  
  
  //Set Name of ROOTfile containing polygon Cuts
  TString outCutFile;
  outCutFile=Form("cuts/XpFpXFp_%s_%d_cut.root",OpticsID.Data(),FileID);
  //------------------ Flag 1 ----> Set the Polygon Cut ----------------------

  
  
  TCanvas *histView_Cut; 
  histView_Cut= new TCanvas("histView_Cut","cut",900,500);
  histView_Cut->Divide(1,1);
  
  //Loop over focal plane variables
  histView_Cut->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogz();    
  
	
  cout << "outCutFile =  " << outCutFile << endl;
  TFile fcut(outCutFile,"UPDATE");
  
  TCutG*t;
  Double_t xtar=0,xptar=0;
  vector<vector<vector<Double_t> > > vxfp;
  vector<vector<vector<Double_t> > > vxpfp;
  vector<vector<vector<TLine*> > > vline;
  vxfp.resize(nxsieve);
  vxpfp.resize(nxsieve);
  vline.resize(nxsieve);
  for  (Int_t nc=0;nc<nxsieve;nc++) {
    vxfp[nc].resize(ndelcut);
    vxpfp[nc].resize(ndelcut);
    vline[nc].resize(ndelcut);
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      vxfp[nc][nd].resize(NumFoil);
      vxpfp[nc][nd].resize(NumFoil);
      vline[nc][nd].resize(NumFoil);
    }
  }
  /*
  Double_t zdis_sieve = 1.172;//to sieve front
  Double_t ytar,yptar;
  for  (Int_t nyscol=0;nyscol<13;nyscol++) {
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      for  (Int_t nf=0;nf<NumFoil;nf++) {
	Double_t DelCent = delcut[nd];//(delcut[nd+1]+delcut[nd])/2;
	yptar = (ys_cent[0]-ztar_foil[nf]*TMath::Sin(CentAngle))/(zdis_sieve-ztar_foil[nf]*TMath::Cos(CentAngle));
	xptar = (xs_cent[nyscol])/(zdis_sieve-ztar_foil[nf]*TMath::Cos(CentAngle));
	ytar = -ztar_foil[nf]*(-TMath::Sin(CentAngle)+yptar*TMath::Cos(CentAngle));
	
	Double_t yfp = 0.0,ypfp=0.0,xpfp=0.0,xfp=0.0;
	Double_t yfplo = 0.0,ypfplo=0.0;
	Double_t yfphi = 0.0,ypfphi=0.0;
	Double_t xfplo = 0.0,xpfplo=0.0;
	Double_t xfphi = 0.0,xpfphi=0.0;
	Double_t etemp;
	for( int icoeffold=0; icoeffold<num_recon_terms; icoeffold++ ){
	  etemp= 
	    pow( xtar, xtarexpon[icoeffold] ) * 
	    pow( xptar*1000, xptarexpon[icoeffold] ) * 
	    pow( ytar, ytarexpon[icoeffold] ) * 
	    pow( yptar*1000, yptarexpon[icoeffold] ) * 
	    pow( DelCent, deltaexpon[icoeffold] );
	  xfp += xfpcoeffs[icoeffold] * etemp;
	  xpfp += xpfpcoeffs[icoeffold] * etemp;
	  yfp += yfpcoeffs[icoeffold] * etemp;
	  ypfp += ypfpcoeffs[icoeffold] * etemp;
	  etemp= 
	    pow( xtar, xtarexpon[icoeffold] ) * 
	    pow( xptar*1000, xptarexpon[icoeffold] ) * 
	    pow( ytar, ytarexpon[icoeffold] ) * 
	    pow( yptar*1000, yptarexpon[icoeffold] ) * 
	    pow(delcut[nd] , deltaexpon[icoeffold] );
	  xfplo += xfpcoeffs[icoeffold] * etemp;
	  xpfplo += xpfpcoeffs[icoeffold] *etemp; 
	  etemp= 
	    pow( xtar, xtarexpon[icoeffold] ) * 
	    pow( xptar*1000, xptarexpon[icoeffold] ) * 
	    pow( ytar, ytarexpon[icoeffold] ) * 
	    pow( yptar*1000, yptarexpon[icoeffold] ) * 
	    pow(delcut[nd+1] , deltaexpon[icoeffold] );
	  xfphi += xfpcoeffs[icoeffold] * etemp;
	  xpfphi += xpfpcoeffs[icoeffold] *etemp; 
	} // for icoeffold loop
	vxfp[nyscol][nd][nf]= xfp;
	vxpfp[nyscol][nd][nf] = xpfp;
	vline[nyscol][nd][nf]= new TLine(xpfplo/1000.,xfplo,xpfphi/1000.,xfphi);
	vline[nyscol][nd][nf]->SetLineColor(1);
	vline[nyscol][nd][nf]->SetLineWidth(3);
	// cout << " full matrix xfp = " << vxfp[nyscol][nd][nf]<< " xpfp = " << vxpfp[nyscol][nd][nf]/1000. << endl;
	//
      }}}
*/
  
  //
  Int_t yscol=0;
  yscol=0;
  Int_t nd=0;
  Int_t nf=0;
  Int_t nloop=0;
  while( nd!=-1 ) {
    cout << " Which ndelta region (nd=-1 to quit) ? " << " nd = " << nd << endl;
    cin >> nd;
    if (nd >=ndelcut) nd=0;
    if (nd !=-1)cout << " Which foil ? " << " nf = " << nf<< endl;
    if (nd !=-1) cin >> nf;
    if (nf >=NumFoil)nf=0;
    Double_t DelCent = delcut[nd];//(delcut[nd+1]+delcut[nd])/2;
    nloop=-1;
    if (nd !=-1) nloop=0;
    while (nloop !=-1 ) {
      fhistroot->cd();
      histView_Cut->Clear();
      hXpFpXFp_DelCut[nf][nd]->Draw("colz");
      hXpFpXFp_DelCut[nf][nd]->SetMinimum(hist_minZ);	
      hXpFpXFp_DelCut[nf][nd]->GetYaxis()->SetRangeUser(AxisRange[nd][0],AxisRange[nd][1]);
      hXpFpXFp_DelCut[nf][nd]->GetXaxis()->SetRangeUser(AxisRange[nd][2],AxisRange[nd][3]);
      histView_Cut->Update();
      for  (Int_t nys=0;nys<nxsieve;nys++) {
	//vline[nys][nd][nf]->Draw("same");
	//TText* ystext = new TText(vxpfp[nys][nd][nf]/1000.,AxisRange[nd][1],Form("P%d",nys));
	//ystext->SetTextColor(2);
	//ystext->Draw();
      }
      histView_Cut->Update();
      fcut.cd();
      for  (Int_t nys=0;nys<nxsieve;nys++) {
	hname_cut= Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d;1",nys,nf,nd);
	t=(TCutG*)gROOT->FindObject(hname_cut);
	//       if(!t){cout << "No cut  =  " <<hname_cut  << endl;}
 	if(t) {
	  //cout << " draw cut = " <<hname_cut<< endl;
	  fcut.cd();
	  t->Draw("same");
	  t->SetLineColor(1);
	  //t->Print();
	  Double_t xcut,ycut;
	  t->GetPoint(0,xcut,ycut);
	  TText* ystext = new TText(xcut,ycut,Form("%d",nys));
	  ystext->Draw();
	  histView_Cut->Update();
	}
      }
      cout <<" Action ( 0 ( set cut) , -1 (pick next nd and foil), -10 delete cut) "  << endl;
      cin >> nloop ;
      if (nloop == -100) return;
      if (!(nloop == -10 || nloop==0 || nloop==-1 || nloop==-100)) continue;//return;
      if (nloop==-10 || nloop ==0) {
	cout << " Which xsieve hole ? " << " yscol = " << yscol << endl;
	cin >> yscol;
	if (yscol >=nxsieve)yscol=0;
	hname_cut= Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d;1",yscol,nf,nd);
      }
      if(nloop==-10) {
	fcut.cd();
	// fcut.ls();
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
	  hname_cut= Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",yscol,nf,nd);
	  t=(TCutG*)gROOT->FindObject(hname_cut);
	  gROOT->Remove(t);
	}
	
	Double_t minx=0;
	Double_t miny=0;
	Double_t maxx=0;
	Double_t maxy=0;
	Int_t ans=-1;
	cout << " Change the axis? ( -1=NO)"<< endl;
	cin >> ans;
	if (ans!=-1) {
	  cout << " What is the X min ? "<< endl;
	  cin >> minx;
	  cout << " What is the X max ? "<< endl;
	  cin >> maxx;
	  cout << " What is the Y min ? "<< endl;
	  cin >> miny;
	  cout << " What is the Y max ? "<< endl;
	  cin >> maxy;
	  fhistroot->cd();
	  histView_Cut->Clear();
	  hXpFpXFp_DelCut[nf][nd]->Draw("colz");
	  hXpFpXFp_DelCut[nf][nd]->SetMinimum(hist_minZ);	
	  hXpFpXFp_DelCut[nf][nd]->GetYaxis()->SetRangeUser(miny,maxy);
	  hXpFpXFp_DelCut[nf][nd]->GetXaxis()->SetRangeUser(minx,maxx);
	  histView_Cut->Update();
	  fcut.cd();
	  for  (Int_t nys=0;nys<nxsieve;nys++) {
	    hname_cut= Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d;1",nys,nf,nd);
	    t=(TCutG*)gROOT->FindObject(hname_cut);
	    //       if(!t){cout << "No cut  =  " <<hname_cut  << endl;}
	    if(t) {
	      //cout << " draw cut = " <<hname_cut<< endl;
	      fcut.cd();
	      t->Draw("same");
	      t->SetLineColor(1);
	      //t->Print();
	      Double_t xcut,ycut;
	      t->GetPoint(0,xcut,ycut);
	      TText* ystext = new TText(xcut,ycut,Form("%d",nys));
	      ystext->Draw();
	      histView_Cut->Update();
	    }
	  }
	}
	TCutG*cutg=(TCutG*)gPad->WaitPrimitive("CUTG","CutG");
	histView_Cut->Update();
	TCutG*tmpg= (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
	hname_cut= Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",yscol,nf,nd);
	TCutG*mycutg=(TCutG*)(tmpg->Clone(hname_cut));
	fcut.cd();
	mycutg->Write("",TObject::kOverwrite);
	//mycutg->Print();
	// mycutg->Draw();
	// histView_Cut->Update();
      }
    }
  }
}


