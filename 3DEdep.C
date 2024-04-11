{
   
    TCanvas* c1 = new TCanvas ("c1","",20,20,800,600);
    gStyle->SetPalette(kThermometer);//;kBlackBody;kThermometer;kMint
   
  // Draw histos filled by Geant4 simulation 
  //   
    TFile* f = new TFile("BG_t0.root");

    TNtuple* ntuple;
    ntuple = (TNtuple*)f->Get("test"); 

    ntuple->Draw("xposition:yposition:edep","colz");
    // TH1F *hist1 = (TH1F*)gDirectory->Get("hist1");
    // hist1->Fit("gaus");

}
