#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLegend.h>
#include <TTree.h>

#include <iostream>

using namespace std;
#include "CommandLine.h"
#include "Messenger.h"
#include "ProgressBar.h"
#include "include/MITHIG_CMSStyle.h"

void divideByWidth(TH1* hist) {
  if (!hist) {
    std::cerr << "Error: Null histogram pointer passed to divideByWidth function." << std::endl;
    return;
  }

  
  int nBins = hist->GetNbinsX();
  for (int i = 1; i <= nBins; ++i) {
  
    double binContent = hist->GetBinContent(i);
    double binError = hist->GetBinError(i);
    double binWidth = hist->GetBinWidth(i);
  

    if (binWidth != 0) {
  
      hist->SetBinContent(i, binContent / binWidth);
      hist->SetBinError(i, binError / binWidth);
    } else {
      std::cerr << "Warning: Bin width is zero for bin " << i << ". Skipping division for this bin." << std::endl;
    }
  }
}

//============================================================//
// Plotting
//============================================================//
int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);
  std::vector<string> inputFilenames = CL.GetStringVector("Input", "input.root");    // Input file
  std::vector<string> SystLabels = CL.GetStringVector("SystLabels", "HFAND15");    // Input file
  std::vector<int> ZBPrescale = CL.GetIntVector("ZBPrescale", "");    // Input file
  string output = CL.Get("Output", "output"); // Output file
  bool doRatio = CL.GetBool("DoRatio", false);

  std::vector<Int_t> colors = {cmsBlue, cmsRed, cmsPurple, cmsTeal, cmsYellow, cmsGray};

  if (inputFilenames.size() != SystLabels.size()) {
    cerr << "Error: The number of input files and labels must match." << endl;
    return 1;
  }

  if (ZBPrescale.size() != inputFilenames.size()) {
    cerr << "Error: The number of ZB prescale flags must match the number of input files." << endl;
    return 1;
  }

  if (inputFilenames.size() > colors.size()) {
    cerr << "Error: More input files than available colors." << endl;
    return 1;
  }

  std::vector<TH1D*> hTrkPt;
  std::vector<TH1D*> hTrkPtRatio;
  std::vector<TH1D*> hTrkEta;
  std::vector<TH1D*> hTrkEtaRatio;
  std::vector<TH1D*> hZBPrescale;

  TCanvas *cTrkPt = new TCanvas("cTrkPt", "Charged Hadron RAA", 800, 600);
  TCanvas *cTrkEta = new TCanvas("cTrkEta", "Charged Hadron RAA", 800, 600);
  TLegend* legendSyst = new TLegend(0.55, 0.67, 0.85, 0.87);
  legendSyst->SetBorderSize(0);
  legendSyst->SetFillStyle(0);
  legendSyst->SetTextFont(42);
  legendSyst->SetTextSize(0.035);

  std::vector<TFile*> inputFiles;

  for (int iFile = 0; iFile < inputFilenames.size(); iFile++) {
  
    inputFiles.push_back(new TFile(inputFilenames[iFile].c_str(), "READ"));
  
    if (!inputFiles[iFile]->IsOpen()) {
      cerr << "Error: Could not open input file " << inputFilenames[iFile] << endl;
      return 1;
    }
  
    hTrkEta.push_back((TH1D*)inputFiles[iFile]->Get("hTrkEta")->Clone(Form("hTrkEta_%i", (iFile))));
    hTrkPt.push_back((TH1D*)inputFiles[iFile]->Get("hTrkPt")->Clone(Form("hTrkPt_%i", (iFile))));
    hZBPrescale.push_back((TH1D*)inputFiles[iFile]->Get("hZBPrescale")->Clone(Form("hZBPrescale_%i", (iFile))));
  }

  cTrkPt->cd();
  if (doRatio) {
    cTrkPt->Divide(1, 2);
    cTrkPt->cd(1);
  } 
  // Get the canvas pad to pass to other functions
  TPad* padTrkPt;
  if (doRatio) {
    padTrkPt = (TPad*) cTrkPt->GetPad(1);
  } else {
    padTrkPt = (TPad*) cTrkPt->GetPad(0);
  }
  padTrkPt->cd();
  padTrkPt->SetLogy();
  padTrkPt->SetLogx();
  padTrkPt->SetTicks();
  if (doRatio) {
    padTrkPt->SetPad(0,0.3,1,1);
    padTrkPt->SetLeftMargin(0.15);
    padTrkPt->SetRightMargin(0.05);
    padTrkPt->SetBottomMargin(0);
    padTrkPt->SetTopMargin(0.1);
  }

  // >>> Apply the CMS TDR style <<<
  SetTDRStyle();

  for (int iFile = 0; iFile < hTrkPt.size(); iFile++) {
    divideByWidth(hTrkPt[iFile]);
    if (ZBPrescale[iFile] > 0) {
      hTrkPt[iFile]->Scale(hZBPrescale[iFile]->GetBinContent(1) / hZBPrescale[iFile]->GetBinContent(2));
    }
    hTrkPt[iFile]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hTrkPt[iFile]->GetYaxis()->SetTitle("dN/dp_{T}");
    hTrkPt[iFile]->GetXaxis()->SetRangeUser(1, 50);
    if(doRatio) {
      hTrkPt[iFile]->GetXaxis()->SetLabelSize(0.10*0.3/0.7);
      hTrkPt[iFile]->GetXaxis()->SetTitleSize(0.12*0.3/0.7);
      hTrkPt[iFile]->GetYaxis()->SetTitleOffset(0.3*7/3);
      hTrkPt[iFile]->GetYaxis()->SetLabelSize(0.10*0.3/0.7);
      hTrkPt[iFile]->GetYaxis()->SetTitleSize(0.12*0.3/0.7);
      hTrkPt[iFile]->GetYaxis()->ChangeLabel(1, -1, 0);
    }
    hTrkPt[iFile]->SetLineColor(colors[iFile]);
    hTrkPt[iFile]->SetLineWidth(3);
    if (iFile>0) hTrkPt[iFile]->Draw("SAME");
    else hTrkPt[iFile]->Draw();
    
    legendSyst->AddEntry(hTrkPt[iFile], SystLabels[iFile].c_str(), "l");
  }
  legendSyst->Draw();

  AddCMSHeader(
    padTrkPt,
    "Internal",
    false
  );
  AddUPCHeader(padTrkPt, "9.6 TeV", "Run 394153 OO");
  padTrkPt->Update();

  if (doRatio) {
    cTrkPt->cd(2);
    gPad->SetPad(0,0,1,0.3);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.28);
    gPad->SetTopMargin(0);

    gPad->SetLogx();
    //gPad->SetLogy();
    gPad->SetTicks();
    // Draw the ratio
    for (int iFile = 0; iFile < hTrkPt.size(); iFile++) {
      if (iFile > 0) {
        hTrkPtRatio.push_back((TH1D*)hTrkPt[iFile]->Clone(Form("hTrkPtRatio_%i", (iFile))));
        hTrkPtRatio[iFile]->GetYaxis()->SetTitle("Ratio");
        hTrkPtRatio[iFile]->Divide(hTrkPt[0]);
        hTrkPtRatio[iFile]->GetYaxis()->SetRangeUser(0.93, 1.07);
        hTrkPtRatio[iFile]->GetYaxis()->ChangeLabel(-1, -1, 0);
        hTrkPtRatio[iFile]->GetXaxis()->SetLabelSize(0.10);
        hTrkPtRatio[iFile]->GetXaxis()->SetTitleSize(0.12);
        hTrkPtRatio[iFile]->GetYaxis()->SetLabelSize(0.10);
        hTrkPtRatio[iFile]->GetYaxis()->SetTitleSize(0.12);
        hTrkPtRatio[iFile]->GetYaxis()->SetTitleOffset(0.3);
        hTrkPtRatio[iFile]->GetXaxis()->SetTitleOffset(0.9);
        //hTrkPtRatio[iFile]->GetYaxis()->SetNdivisions(505);
        if (iFile == 1) {
          hTrkPtRatio[iFile]->Draw();
        } else {
          hTrkPtRatio[iFile]->Draw("SAME");
        }
      } else {
        hTrkPtRatio.push_back((TH1D*)hTrkPt[iFile]->Clone(Form("hTrkPtRatio_%i", (iFile))));
      }
    }
  }

  cTrkPt->SaveAs(Form("trkPtSystOverlay_%s.pdf", output.c_str()));

  cTrkEta->cd();
  if (doRatio) {
    cTrkEta->Divide(1, 2);
    cTrkEta->cd(1);
  }
  // Get the canvas pad to pass to other functions
  TPad* padTrkEta;
  if (doRatio) {
    padTrkEta = (TPad*) cTrkEta->GetPad(1);
  } else {
    padTrkEta = (TPad*) cTrkEta->GetPad(0);
  }
  padTrkEta->cd();
  padTrkEta->SetLogy();
  if (doRatio) {
    padTrkEta->SetPad(0,0.3,1,1);
    padTrkEta->SetLeftMargin(0.15);
    padTrkEta->SetRightMargin(0.05);
    padTrkEta->SetBottomMargin(0);
    padTrkEta->SetTopMargin(0.1);
  }

  // >>> Apply the CMS TDR style <<<
  SetTDRStyle();

  for (int iFile = 0; iFile < hTrkEta.size(); iFile++) {
    divideByWidth(hTrkEta[iFile]);
    if (ZBPrescale[iFile] > 0) {
      hTrkEta[iFile]->Scale(hZBPrescale[iFile]->GetBinContent(1) / hZBPrescale[iFile]->GetBinContent(2));
    }
    hTrkEta[iFile]->GetXaxis()->SetTitle("#eta");
    hTrkEta[iFile]->GetYaxis()->SetTitle("dN/d#eta");
    hTrkEta[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
    if(doRatio) {
      hTrkEta[iFile]->GetXaxis()->SetLabelSize(0.10*0.3/0.7);
      hTrkEta[iFile]->GetXaxis()->SetTitleSize(0.12*0.3/0.7);
      hTrkEta[iFile]->GetYaxis()->SetTitleOffset(0.3*7/3);
      hTrkEta[iFile]->GetYaxis()->SetLabelSize(0.10*0.3/0.7);
      hTrkEta[iFile]->GetYaxis()->SetTitleSize(0.12*0.3/0.7);
      //hTrkEta[iFile]->GetYaxis()->ChangeLabel(1, -1, 0);
    }
    hTrkEta[iFile]->SetLineColor(colors[iFile]);
    hTrkEta[iFile]->SetLineWidth(3);
    if (iFile>0) hTrkEta[iFile]->Draw("SAME");
    else hTrkEta[iFile]->Draw();
  }
  legendSyst->Draw();

  AddCMSHeader(
    padTrkEta,
    "Internal",
    false
  );
  AddUPCHeader(padTrkEta, "9.6 TeV", "Run 394153 OO");
  padTrkEta->Update();
  if (doRatio) {
    cTrkEta->cd(2);
    gPad->SetPad(0,0,1,0.3);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.28);
    gPad->SetTopMargin(0);

    //gPad->SetLogx();
    //gPad->SetLogy();
    gPad->SetTicks();
    // Draw the ratio
    for (int iFile = 0; iFile < hTrkEta.size(); iFile++) {
      if (iFile > 0) {
        hTrkEtaRatio.push_back((TH1D*)hTrkEta[iFile]->Clone(Form("hTrkEtaRatio_%i", (iFile))));
        hTrkEtaRatio[iFile]->GetYaxis()->SetTitle("Ratio");
        hTrkEtaRatio[iFile]->Divide(hTrkEta[0]);
        hTrkEtaRatio[iFile]->GetYaxis()->SetRangeUser(0.93, 1.07);
        hTrkEtaRatio[iFile]->GetYaxis()->ChangeLabel(-1, -1, 0);
        hTrkEtaRatio[iFile]->GetXaxis()->SetLabelSize(0.10);
        hTrkEtaRatio[iFile]->GetXaxis()->SetTitleSize(0.12);
        hTrkEtaRatio[iFile]->GetYaxis()->SetLabelSize(0.10);
        hTrkEtaRatio[iFile]->GetYaxis()->SetTitleSize(0.12);
        hTrkEtaRatio[iFile]->GetYaxis()->SetTitleOffset(0.3);
        hTrkEtaRatio[iFile]->GetXaxis()->SetTitleOffset(0.9);
        //hTrkEtaRatio[iFile]->GetYaxis()->SetNdivisions(505);
        if (iFile == 1) {
          hTrkEtaRatio[iFile]->Draw();
        } else {
          hTrkEtaRatio[iFile]->Draw("SAME");
        }
      } else {
        hTrkEtaRatio.push_back((TH1D*)hTrkEta[iFile]->Clone(Form("hTrkEtaRatio_%i", (iFile))));
      }
    }
  }
  cTrkEta->SaveAs(Form("trkEtaSystOverlay_%s.pdf", output.c_str()));

  for (int iFile = 0; iFile < inputFiles.size(); iFile++) {
    inputFiles[iFile]->Close();
  }
  return 0;
}