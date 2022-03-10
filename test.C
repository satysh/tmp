Double_t muI(Double_t *x, Double_t *par) {

  Double_t p=x[0]; // p

  Double_t ksi=p*x[1]; // p*cost

  const Double_t c1=0.00253;
  const Double_t c2=0.2455;
  const Double_t c3=1.288;
  const Double_t c4=-0.2555;
  const Double_t c5=0.0209;

  Double_t lksi=TMath::Log10(ksi);

  Double_t Iv = par[0]*c1*pow(ksi,-1.*(c2+c3*lksi+c4*lksi*lksi+c5*lksi*lksi*lksi));

  return Iv*x[1]*x[1]*x[1];

}

void test()
{
	Int_t nEvents = 10000;
  gRandom->SetSeed(0);
  
  const Int_t npar = 1;
  
  TF2 *f2 = new TF2("f2",muI,0,2000,0,1, npar);
  Double_t f2params[npar] = {0.8915};
  f2->SetParameters(f2params);
  f2->SetNpx(5000);
  //f2->Draw();
  
  Double_t cost, p;
  Double_t theta;
  Double_t phi;
  Float_t rnd;
  TTree *tree = new TTree("tree", "tree");
  tree->Branch("p", &p);
  tree->Branch("cost", &cost);
  tree->Branch("theta", &theta);
  tree->Branch("phi", &phi);
  tree->Branch("rnd", &rnd);
  for (Int_t i=0; i<nEvents; i++) {
    f2->GetRandom2(p, cost );
    theta = TMath::RadToDeg()*acos(cost);
    rnd = gRandom->Rndm();
    phi = TMath::RadToDeg()*2.*TMath::Pi()*rnd;
    tree->Fill();
  }
  tree->StartViewer();

/*
  TCanvas *canv = new TCanvas("canv", "canv");
  canv->Divide(1, 2);
  canv->cd(1);
  tree->Draw("p");
  canv->cd(2);
  tree->Draw("cost");
*/  
}