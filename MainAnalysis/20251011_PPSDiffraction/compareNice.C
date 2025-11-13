#include "MITHIG_CMSStyle.h"
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>

// Compare histograms from two different output files
void draw(TString tag1 = "ppsTag_0_zdcPlusVeto_0_HFPlusGapVeto_1",
          TString tag2 = "ppsTag_0_zdcPlusVeto_0_HFPlusGapVeto_1",
          TString legendName1 = "No cuts",
          TString legendName2 = "With cuts",
          TString histname = "hNtrk",
          float xmin = 0, float xmax = 100,
          float ymin = 1e-6, float ymax = 100,  // ymin > 0 for logY
          bool logy = true, bool logx = false) {

  TFile *f1 = TFile::Open(Form("output_loop_%s.root", tag1.Data()));
  TFile *f2 = TFile::Open(Form("output_loop_%s.root", tag2.Data()));
  TH1F *h1 = (TH1F*)f1->Get(histname);
  TH1F *h2 = (TH1F*)f2->Get(histname);

  TCanvas *c1 = new TCanvas(Form("c1_%s", histname.Data()), "Comparison", 600, 600);
  TPad* pad1 = (TPad*) c1->GetPad(0);
  pad1->cd();
  SetTDRStyle();

  if (logy) pad1->SetLogy();
  if (logx) pad1->SetLogx();

  // Normalize safely
  double i1 = h1->Integral();
  double i2 = h2->Integral();
  if (i1 > 0) h1->Scale(1.0/i1);
  if (i2 > 0) h2->Scale(1.0/i2);

  // X range zoom is fine with SetRangeUser on X
  h1->GetXaxis()->SetRangeUser(xmin, xmax);

  // Y range must use SetMinimum/SetMaximum on the FIRST drawn hist
  if (logy && ymin <= 0) ymin = 1e-6; // enforce >0 on log pad
  h1->SetMinimum(ymin);
  h1->SetMaximum(ymax);

  SetTH1Fonts(h1, 1.2);
  SetTH1Fonts(h2, 1.2);
  h1->SetLineColor(cmsRed);
  h1->SetLineWidth(2);
  h2->SetLineColor(cmsBlue);
  h2->SetLineWidth(2);

  h1->Draw("HIST");         // first draw defines the frame
  h2->Draw("HIST SAME");    // no autoscale now

  TLegend *legend = new TLegend(0.2, 0.75, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(h1, legendName1, "l");
  legend->AddEntry(h2, legendName2, "l");
  legend->Draw();

  AddCMSHeader(pad1, "Internal", false);
  AddUPCHeader(pad1, "9.62 TeV", "pO");
  pad1->Update();

  c1->Modified(); c1->Update();
  c1->SaveAs(Form("plots/comparison_%s_vs_%s_%s.pdf",
                  tag1.Data(), tag2.Data(), histname.Data()));
}

void compare(TString string1 = "ppsTag_0_zdcPlusVeto_0_HFPlusGapVeto_0",
             TString string2 = "ppsTag_0_zdcPlusVeto_1_HFPlusGapVeto_1",
             TString legendName1 = "No cuts",
             TString legendName2 = "With cuts") {
  draw(string1, string2, legendName1, legendName2, "hNtr", 0, 100, 1e-4, 2, true, false);
  draw(string1, string2, legendName1, legendName2, "hEtaCharged", -3, 3, 1e-4, 0.5, true, false);
  draw(string1, string2, legendName1, legendName2, "hZDCm", 0, 6000, 1e-5, 2, true, false);
  draw(string1, string2, legendName1, legendName2, "hZDCp", 0, 6000, 1e-5, 2, true, false);
}

int main() {
  compare("ppsTag_0_zdcPlusVeto_0_HFPlusGapVeto_0",
          "ppsTag_0_zdcPlusVeto_1_HFPlusGapVeto_1",
          "No cuts",
          "With ZDC+ veto and HF+ gap veto");
  return 0;
}