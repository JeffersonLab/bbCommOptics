void ratePlotter(){
  TString inputroot = "hist/Optics_gen_c1_sieve4_-1_fit_tree.root";
  TString outputpdf= "output/sieve_rates_c1_s4";
  
  //setup foils
  const int nfoil = 10;
  vector<double> zfoil = {0.2,0.15,0.1,0.05,0.0,-0.05,-0.1,-0.15,-0.2,-0.25};

  //setup sieve
  const int nxsieve = 19;//17;13;
  const int nysieve = 13;//11;7;

  vector <Double_t> ys_cent;
  for (Int_t nys=0;nys<nysieve;nys++) {
    //Double_t pos=nys*0.0381-0.0381*3;//old sieve
    //Double_t pos = nys*0.0254-0.0254*5;
    Double_t pos = nys*(7./8.)*0.0254-(7./8.)*0.0254*6;//newest sieve
    cout<<"nys: "<<nys<<" pos: "<<pos<<endl;
    ys_cent.push_back(pos);
  }
  
  vector <Double_t> xs_cent;
  for (Int_t nys=0;nys<nxsieve;nys++) {
    //Double_t pos=xs_cent[nys];//old sieve
    //Double_t pos = nys*(1.+9./16.)*0.0254-(1.+9./16.)*0.0254*8.;
    Double_t pos = nys*(1.25)*0.0254-(1.25)*0.0254*9.;//newest sieve
    xs_cent.push_back(pos);
  } 
   
  /*
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
  */

  //read in the input file
  TFile *fsimc = new TFile(inputroot); 
  TTree *FitTree = (TTree*) fsimc->Get("TFit");
  Double_t  ys,xtar,ztar,xptar,yptar,ytar,delta,xptarT,yptarT,ytarT,ztarT;
  Double_t ysieveT,ysieve, xsieve, xsieveT,weight;
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
  FitTree->SetBranchAddress("weight",&weight);

  //make rate vector to count
  vector<vector<vector<double>>> holeRate;//[foil][ys][xs]
  holeRate.resize(nfoil);
  for (int iz=0; iz<nfoil; iz++){
    holeRate[iz].resize(nysieve);
    for (int iy=0; iy<nysieve; iy++){
      holeRate[iz][iy].resize(nxsieve);
    }
  }

  //loop events and count rate
  Long64_t nentries = FitTree->GetEntries();
  for (int i = 0; i < nentries; i++) {
  //for (int i = 0; i < 1000; i++) {
    FitTree->GetEntry(i);
    int izF=-1;
    int iyF=-1;
    int ixF=-1;
    //cout<<"weight: "<<weight<<" ztar: "<<ztarT<<" ysieveT: "<<ysieveT<<" xsieveT: "<<xsieveT<<endl;
    for (int iz=0; iz<nfoil; iz++){
      if (abs(ztarT-zfoil[iz])<0.01) izF=iz;
      for (int iy=0; iy<nysieve; iy++){
	if (abs(ysieveT-ys_cent[iy])<0.005) iyF=iy;
	for (int ix=0; ix<nxsieve; ix++){
	  if (abs(xsieveT-xs_cent[ix])<0.005) ixF=ix;
	}
      }
    }
    if (izF!=-1 && iyF!=-1 && ixF!=-1){
      holeRate[izF][iyF][ixF] += weight;
    }
  }

  //plot the rates per foil per hole
  TH2F *h_rateAllZ = new TH2F("h_rateAllZ","rate for all foils combined; y sieve; x sieve",nysieve,0,nysieve,nxsieve,0,nxsieve);

  vector <TH2F*> h_rates;
  h_rates.resize(nfoil);

  for (int iz=0; iz<nfoil; iz++){
    h_rates[iz] = new TH2F(Form("h_rates_%d",iz),Form("rates, Foil at %2.2f;y sieve; x sieve",zfoil[iz]),nysieve,0,nysieve,nxsieve,0,nxsieve);
    for (int iy=0; iy<nysieve; iy++){
      for (int ix=0; ix<nxsieve; ix++){
	Int_t bin = h_rates[iz]->GetBin(iy+1,ix+1,0);
	h_rates[iz]->SetBinContent(bin,holeRate[iz][iy][ix]);
	h_rateAllZ->AddBinContent(bin,holeRate[iz][iy][ix]);
      }
    }
    double rPerFoil = h_rates[iz]->Integral();
    cout<<"rate for foil at "<<zfoil[iz]<<": "<<rPerFoil<<endl;
  }

  double totalRate = h_rateAllZ->Integral();
  cout<<"total rate: "<<totalRate<<endl;

  
  //setup output pdf
  TCanvas* can2d[nfoil+1];
  gStyle->SetOptStat(0);
  for  (Int_t nc=0;nc<nfoil;nc++) {
    can2d[nc] = new TCanvas(Form("Can2d_%d",nc),Form("foil %3.2f",zfoil[nc]), 700,700);
    //h_rates[nc]->SetMaximum(1.);
    //h_rates[nc]->SetMinimum(0.);
    h_rates[nc]->Draw("colz");
       
    TString end = ".pdf";
    if (nc==0) end=".pdf(";
    //if (nc==nfoil-1) end=".pdf)";
    can2d[nc]->Print(outputpdf+end);
  }
  can2d[nfoil] = new TCanvas(Form("Can2d_%d",nfoil),"all foils", 700,700);
  h_rateAllZ->SetEntries(1);
  h_rateAllZ->Draw("colz");
  can2d[nfoil]->Print(outputpdf+".pdf)");
  
}










  







