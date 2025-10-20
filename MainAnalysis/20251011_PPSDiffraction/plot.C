#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm> // for std::min
#include <iostream>
#include <vector>

void plot(char *filename) {

  TFile *f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open file!\n";
    return;
  }

  TH1F *hNtrNoCond = (TH1F *)f->Get("hNtrNoCond");
  TH1F *hNtrPPS = (TH1F *)f->Get("hNtrPPS");

  TCanvas *cChargedTracks = new TCanvas("cChargedTracks", "cChargedTracks", 1000, 850);
  cChargedTracks->SetLogy();
  hNtrNoCond->Rebin(4);
  hNtrPPS->Rebin(4);
  hNtrNoCond->Scale(1. / hNtrNoCond->Integral());
  hNtrNoCond->SetLineColor(kRed);
  hNtrPPS->Scale(1. / hNtrPPS->Integral());
  hNtrPPS->SetLineColor(kBlue);
  hNtrNoCond->GetXaxis()->SetTitle("Charged hadron multiplicity, #eta|< 2.4");
  hNtrNoCond->GetXaxis()->SetRangeUser(0, 300);
  hNtrPPS->GetXaxis()->SetRangeUser(0, 300);
  hNtrNoCond->SetMaximum(1.0);
  hNtrNoCond->Draw("same");
  hNtrPPS->Draw("same");
  TLegend *legend = new TLegend(0.4, 0.7, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.035);
  legend->AddEntry(hNtrNoCond, "No requirements", "l");
  legend->AddEntry(hNtrPPS, "At least 1 proton in PPS", "l");
  legend->Draw();
  cChargedTracks->SaveAs("cChargedTracks.pdf");

  TH1F *hFSC2bottomM_adcPPS = (TH1F *)f->Get("hFSC2bottomM_adcPPS");
  TH1F *hFSC2bottomM_adcNoCond = (TH1F *)f->Get("hFSC2bottomM_adcNoCond");

  TCanvas *cFSC2bottomM = new TCanvas("cFSC2bottomM", "cFSC2bottomM", 1000, 850);
  cFSC2bottomM->SetLogy();
  hFSC2bottomM_adcPPS->Rebin(4);
  hFSC2bottomM_adcNoCond->Rebin(4);
  hFSC2bottomM_adcPPS->Scale(1. / hFSC2bottomM_adcPPS->Integral());
  hFSC2bottomM_adcPPS->SetLineColor(kBlue);
  hFSC2bottomM_adcNoCond->Scale(1. / hFSC2bottomM_adcNoCond->Integral());
  hFSC2bottomM_adcNoCond->SetLineColor(kRed);
  hFSC2bottomM_adcPPS->GetXaxis()->SetTitle("FSC2 bottom ADC");
  hFSC2bottomM_adcPPS->GetXaxis()->SetRangeUser(0, 2000);
  hFSC2bottomM_adcNoCond->GetXaxis()->SetRangeUser(0, 2000);
  hFSC2bottomM_adcPPS->SetMaximum(1.0);
  hFSC2bottomM_adcPPS->Draw("same");
  hFSC2bottomM_adcNoCond->Draw("same");
  TLegend *legend2 = new TLegend(0.4, 0.7, 0.9, 0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry(hFSC2bottomM_adcPPS, "At least 1 proton in PPS", "l");
  legend2->AddEntry(hFSC2bottomM_adcNoCond, "No requirements", "l");
  legend2->Draw();
  cFSC2bottomM->SaveAs("cFSC2bottomM.pdf");

  TH1F *hFSC2bottomM_chargefCNoCond = (TH1F *)f->Get("hFSC2bottomM_chargefCNoCond");
  TH1F *hFSC2bottomM_chargefCPPS = (TH1F *)f->Get("hFSC2bottomM_chargefCPPS");

  TCanvas *cFSC2bottomM_chargefC = new TCanvas("cFSC2bottomM_chargefC", "cFSC2bottomM_chargefC", 1000, 850);
  cFSC2bottomM_chargefC->SetLogy();
  hFSC2bottomM_chargefCPPS->GetXaxis()->SetRangeUser(0, 400000);
  hFSC2bottomM_chargefCNoCond->GetXaxis()->SetRangeUser(0, 400000);
  hFSC2bottomM_chargefCPPS->Rebin(10);
  hFSC2bottomM_chargefCNoCond->Rebin(10);
  hFSC2bottomM_chargefCPPS->Scale(1. / hFSC2bottomM_chargefCPPS->Integral());
  hFSC2bottomM_chargefCPPS->SetLineColor(kBlue);
  hFSC2bottomM_chargefCNoCond->Scale(1. / hFSC2bottomM_chargefCNoCond->Integral());
  hFSC2bottomM_chargefCNoCond->SetLineColor(kRed);
  hFSC2bottomM_chargefCPPS->GetXaxis()->SetTitle("FSC2 bottom charge [fC]");

  hFSC2bottomM_chargefCPPS->SetMaximum(4.0);
  hFSC2bottomM_chargefCPPS->Draw();
  hFSC2bottomM_chargefCNoCond->Draw("same");
  TLegend *legend3 = new TLegend(0.4, 0.7, 0.9, 0.9);
  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->SetTextSize(0.035);
  legend3->AddEntry(hFSC2bottomM_chargefCPPS, "At least 1 proton in PPS", "l");
  legend3->AddEntry(hFSC2bottomM_chargefCNoCond, "No requirements", "l");
  legend3->Draw();
  cFSC2bottomM_chargefC->SaveAs("cFSC2bottomM_chargefC.pdf");

  TCanvas *cFSC3bottomleftM = new TCanvas("cFSC3bottomleftM", "cFSC3bottomleftM", 1000, 850);
  TH1F *hFSC3bottomleftM_adcPPS = (TH1F *)f->Get("hFSC3bottomleftM_adcPPS");
  TH1F *hFSC3bottomleftM_adcNoCond = (TH1F *)f->Get("hFSC3bottomleftM_adcNoCond");
  cFSC3bottomleftM->SetLogy();
  hFSC3bottomleftM_adcPPS->Rebin(4);
  hFSC3bottomleftM_adcNoCond->Rebin(4);
  hFSC3bottomleftM_adcPPS->Scale(1. / hFSC3bottomleftM_adcPPS->Integral());
  hFSC3bottomleftM_adcPPS->SetLineColor(kBlue);
  hFSC3bottomleftM_adcNoCond->Scale(1. / hFSC3bottomleftM_adcNoCond->Integral());
  hFSC3bottomleftM_adcNoCond->SetLineColor(kRed);
  hFSC3bottomleftM_adcPPS->GetXaxis()->SetTitle("FSC3 bottom left ADC");
  hFSC3bottomleftM_adcPPS->GetXaxis()->SetRangeUser(0, 2000);
  hFSC3bottomleftM_adcNoCond->GetXaxis()->SetRangeUser(0, 2000);
  hFSC3bottomleftM_adcPPS->SetMaximum(1.0);
  hFSC3bottomleftM_adcPPS->Draw("same");
  hFSC3bottomleftM_adcNoCond->Draw("same");
  TLegend *legend4 = new TLegend(0.4, 0.7, 0.9, 0.9);
  legend4->SetBorderSize(0);
  legend4->SetFillStyle(0);
  legend4->SetTextSize(0.035);
  legend4->AddEntry(hFSC3bottomleftM_adcPPS, "At least 1 proton in PPS", "l");
  legend4->AddEntry(hFSC3bottomleftM_adcNoCond, "No requirements", "l");
  legend4->Draw();
  cFSC3bottomleftM->SaveAs("cFSC3bottomleftM.pdf");

  // FSC3 bottom left charge plots
  TCanvas *cFSC3bottomleftM_chargefC = new TCanvas("cFSC3bottomleftM_chargefC", "cFSC3bottomleftM_chargefC", 1000, 850);
  TH1F *hFSC3bottomleftM_chargefCPPS = (TH1F *)f->Get("hFSC3bottomleftM_chargefCPPS");
  TH1F *hFSC3bottomleftM_chargefCNoCond = (TH1F *)f->Get("hFSC3bottomleftM_chargefCNoCond");
  cFSC3bottomleftM_chargefC->SetLogy();
  hFSC3bottomleftM_chargefCPPS->GetXaxis()->SetRangeUser(0, 400000);
  hFSC3bottomleftM_chargefCNoCond->GetXaxis()->SetRangeUser(0, 400000);
  hFSC3bottomleftM_chargefCPPS->Rebin(10);
  hFSC3bottomleftM_chargefCNoCond->Rebin(10);
  hFSC3bottomleftM_chargefCPPS->Scale(1. / hFSC3bottomleftM_chargefCPPS->Integral());
  hFSC3bottomleftM_chargefCPPS->SetLineColor(kBlue);
  hFSC3bottomleftM_chargefCNoCond->Scale(1. / hFSC3bottomleftM_chargefCNoCond->Integral());
  hFSC3bottomleftM_chargefCNoCond->SetLineColor(kRed);
  hFSC3bottomleftM_chargefCPPS->GetXaxis()->SetTitle("FSC3 bottom left charge [fC]");
  hFSC3bottomleftM_chargefCPPS->SetMaximum(4.0);
  hFSC3bottomleftM_chargefCPPS->Draw();
  hFSC3bottomleftM_chargefCNoCond->Draw("same");
  TLegend *legend5 = new TLegend(0.4, 0.7, 0.9, 0.9);
  legend5->SetBorderSize(0);
  legend5->SetFillStyle(0);
  legend5->SetTextSize(0.035);
  legend5->AddEntry(hFSC3bottomleftM_chargefCPPS, "At least 1 proton in PPS", "l");
  legend5->AddEntry(hFSC3bottomleftM_chargefCNoCond, "No requirements", "l");
  legend5->Draw();
  cFSC3bottomleftM_chargefC->SaveAs("cFSC3bottomleftM_chargefC.pdf");

  // ZDC- plots
  TCanvas *cZDCm = new TCanvas("cZDCm", "cZDCm", 1000, 850);
  TH1F *hZDCmPPS = (TH1F *)f->Get("hZDCmPPS");
  TH1F *hZDCmNoCond = (TH1F *)f->Get("hZDCmNoCond");
  cZDCm->SetLogy();
  hZDCmPPS->Rebin(4);
  hZDCmNoCond->Rebin(4);
  hZDCmPPS->Scale(1. / hZDCmPPS->Integral());
  hZDCmPPS->SetLineColor(kBlue);
  hZDCmNoCond->Scale(1. / hZDCmNoCond->Integral());
  hZDCmNoCond->SetLineColor(kRed);
  hZDCmPPS->GetXaxis()->SetTitle("ZDC-");
  hZDCmPPS->GetXaxis()->SetRangeUser(0, 8000);
  hZDCmNoCond->GetXaxis()->SetRangeUser(0, 8000);
  hZDCmPPS->SetMaximum(1.0);
  hZDCmPPS->Draw();
  hZDCmNoCond->Draw("same");
  TLegend *legend6 = new TLegend(0.4, 0.7, 0.9, 0.9);
  legend6->SetBorderSize(0);
  legend6->SetFillStyle(0);
  legend6->SetTextSize(0.035);
  legend6->AddEntry(hZDCmPPS, "At least 1 proton in PPS", "l");
  legend6->AddEntry(hZDCmNoCond, "No requirements", "l");
  legend6->Draw();
  cZDCm->SaveAs("cZDCm.pdf");

  // ZDC+ plots
  TCanvas *cZDCp = new TCanvas("cZDCp", "cZDCp", 1000, 850);
  TH1F *hZDCpPPS = (TH1F *)f->Get("hZDCpPPS");
  TH1F *hZDCpNoCond = (TH1F *)f->Get("hZDCpNoCond");
  cZDCp->SetLogy();
  hZDCpPPS->Rebin(4);
  hZDCpNoCond->Rebin(4);
  hZDCpPPS->Scale(1. / hZDCpPPS->Integral());
  hZDCpPPS->SetLineColor(kBlue);
  hZDCpNoCond->Scale(1. / hZDCpNoCond->Integral());
  hZDCpNoCond->SetLineColor(kRed);
  hZDCpPPS->GetXaxis()->SetTitle("ZDC+");
  hZDCpPPS->GetXaxis()->SetRangeUser(0, 8000);
  hZDCpNoCond->GetXaxis()->SetRangeUser(0, 8000);
  hZDCpPPS->SetMaximum(1);
  hZDCpPPS->Draw();
  hZDCpNoCond->Draw("same");
  TLegend *legend7 = new TLegend(0.4, 0.7, 0.9, 0.9);
  legend7->SetBorderSize(0);
  legend7->SetFillStyle(0);
  legend7->SetTextSize(0.035);
  legend7->AddEntry(hZDCpPPS, "At least 1 proton in PPS", "l");
  legend7->AddEntry(hZDCpNoCond, "No requirements", "l");
  legend7->Draw();
  cZDCp->SaveAs("cZDCp.pdf");
}

int main(char **argv, char **argc) {
  if (argc < 2) {
    std::cerr << "Please provide the input filename as an argument.\n";
    return 1;
  }
  plot(argv[1]);
  return 0;
}
