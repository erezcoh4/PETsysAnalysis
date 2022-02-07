{
gStyle->SetOptStat("");
TFile *f = new TFile("rossi.root");
TH1F * h1  = (TH1F*)f->Get("1");

Double_t scale = 1/h1->Integral();
h1->Scale(scale);
  
TCanvas *c1 = new TCanvas("c1","c1",800,1000);
c1->cd();
h1->Draw("HIST B");

TH1F * h2 = (TH1F*)f->Get("2");

Double_t scale2 = 1/h2->Integral();
h2->Scale(scale2);
  
TCanvas *c2 = new TCanvas("c1","c1",800,1000);
c2->cd();
h2->Draw("HIST B");
}
