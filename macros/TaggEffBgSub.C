void TaggEffBgSub(TString sBeam, TString sBkg1, TString sBkg2="", Bool_t bFreeScal=false){
  TFile fBeam(sBeam,"READ");
  TFile fBkg1(sBkg1,"READ");  
  gROOT->cd();

  TH1D *hBeamAllHits = (TH1D*)fBeam.Get("TaggerAllHits");
  TH1D *hBeamSingles = (TH1D*)fBeam.Get("TaggerSingles");
  TH1D *hBeamDoubles = (TH1D*)fBeam.Get("TaggerDoubles");
  TH1D *hBeamAccScal = (TH1D*)fBeam.Get("TaggerAccScal");
  TH1D *hBeamLiveTime = (TH1D*)fBeam.Get("LiveTimeScal");

  Double_t dBeamClock = hBeamLiveTime->GetBinContent(hBeamLiveTime->GetXaxis()->FindBin("Clock"));
  Double_t dBeamInhib = hBeamLiveTime->GetBinContent(hBeamLiveTime->GetXaxis()->FindBin("Inhibited"));

  TH1D *hBkg1AccScal = (TH1D*)fBkg1.Get("TaggerAccScal");
  TH1D *hBkg1LiveTime = (TH1D*)fBkg1.Get("LiveTimeScal");

  TH1D *hBackAccScal = (TH1D*)hBkg1AccScal->Clone("hBackAccScal");

  Double_t dBackClock = hBkg1LiveTime->GetBinContent(hBkg1LiveTime->GetXaxis()->FindBin("Clock"));
  Double_t dBackInhib = hBkg1LiveTime->GetBinContent(hBkg1LiveTime->GetXaxis()->FindBin("Inhibited"));

  if(sBkg2 != ""){
    TFile fBkg2(sBkg2,"READ");
    gROOT->cd();

    TH1D *hBkg2AccScal = (TH1D*)fBkg2.Get("TaggerAccScal");
    TH1D *hBkg2LiveTime = (TH1D*)fBkg2.Get("LiveTimeScal");

    hBackAccScal->Add(hBkg2AccScal);
    dBackClock += hBkg2LiveTime->GetBinContent(hBkg2LiveTime->GetXaxis()->FindBin("Clock"));
    dBackInhib += hBkg2LiveTime->GetBinContent(hBkg2LiveTime->GetXaxis()->FindBin("Inhibited"));
  }

  TString sOut = sBeam;
  sOut.ReplaceAll("GoAT","BkgSub");

  TFile fOut(sOut,"RECREATE");

  TH1D *hEffAllHits = (TH1D*)hBeamAllHits->Clone("TaggEffAllHits");
  TH1D *hEffSingles = (TH1D*)hBeamSingles->Clone("TaggEffSingles");
  TH1D *hEffDoubles = (TH1D*)hBeamDoubles->Clone("TaggEffDoubles");
  TH1D *hEffAccScal = (TH1D*)hBeamAccScal->Clone("TaggEffAccScal");
  hEffAccScal->Sumw2();

  if(bFreeScal){
    hEffAccScal->Add(hBackAccScal,(-dBeamClock/dBackClock));
    hEffAccScal->Scale(dBeamInhib/dBeamClock);
  }
  else hEffAccScal->Add(hBackAccScal,(-dBeamInhib/dBackInhib));

  hEffAllHits->Divide(hEffAccScal);
  hEffSingles->Divide(hEffAccScal);
  hEffDoubles->Divide(hEffAccScal);

  ofstream fout;
  TString foutname = sBeam;
  foutname.ReplaceAll("GoAT_", "");
  foutname.ReplaceAll("TaggEff_", "TaggEffUncorr_");
  foutname.ReplaceAll(".root", ".dat");
  fout.open(foutname, ios::ate);
  Int_t nbins = hEffAllHits->GetXaxis()->GetNbins();
  for (int jj=1; jj<=nbins; jj++)
    fout << jj << "\t" << hEffAllHits->GetBinContent(jj) << "\t" << hEffAllHits->GetBinError(jj) << endl;
  fout.close();

  Double_t bins[500], rerr[500]={0.};
  for (int kk=0; kk<nbins; kk++) {
    bins[kk] = kk;
    if (hEffAllHits->GetBinContent(kk)>0)
      rerr[kk] = hEffAllHits->GetBinError(kk)/hEffAllHits->GetBinContent(kk);
  }
  TCanvas *c = new TCanvas("c","c", 800, 500);
  c->cd();
  c->SetGridy();
  TGraph *relerr = new TGraph(nbins, bins, rerr);
  relerr->SetMarkerStyle(20);
  relerr->SetMarkerSize(0.8);
  relerr->GetXaxis()->SetTitle("Channel number");
  relerr->GetYaxis()->SetTitle("Relative error");
  relerr->GetYaxis()->SetRangeUser(-0.05, 0.3);
  relerr->SetName("RelErr");
  relerr->SetTitle("");
  relerr->Draw("AP");
  relerr->Write();
  c->Write();

  Double_t dMean = 0;
  Int_t iMean = 0;
  for(Int_t i=0; i<hEffAllHits->GetNbinsX(); i++){
    if(hEffAllHits->GetBinContent(i+1) > 0){
      dMean += hEffAllHits->GetBinContent(i+1);
      iMean++;
    }
  }

  dMean = 0.015*TMath::Nint(100*dMean/iMean);

  hEffAllHits->GetYaxis()->SetRangeUser(0,dMean);
  hEffSingles->GetYaxis()->SetRangeUser(0,dMean);
  hEffDoubles->GetYaxis()->SetRangeUser(0,dMean);
  
  fOut.Write();
  fOut.Close();

  delete c;
}
