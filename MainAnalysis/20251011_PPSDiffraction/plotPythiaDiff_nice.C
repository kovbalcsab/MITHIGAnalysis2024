#include "MITHIG_CMSStyle.h"
#include "CommandLine.h"
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>

void AddPythiaHeader_custom(
  TPad* pad,
  TString extraText = " + EvtGen",
  bool insideFrame = true, // Set false to put CMS label outside of plot frame
  float textSize = plotTextSize,
  float textOffset = 0.015
) {
  pad->cd();
  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);
  latex.SetTextFont(42);
  latex.SetTextSize(textSize);
  
  if (insideFrame) {
    latex.SetTextAlign(13);
    latex.DrawLatex(
      l + (1 - l - r) * 0.05,
      1 - (t + (1 - t - b) * 0.05) - 2.35 * textSize - textOffset,
      Form("P#scale[0.8]{YTHIA}8%s", extraText.Data())
    );
  }
  else {
    latex.SetTextAlign(11);
    latex.DrawLatex(
      l,
      1 - t + textOffset,
      Form("P#scale[0.8]{YTHIA}8%s", extraText.Data())
    );
  }
}

int main(int argc, char** argv) {
    CommandLine cl(argc, argv);

    std::string inputFile = cl.Get("inputFile", "hardDiffractionUPCs_Op.root");
    TString histName = cl.Get("histName", "hetahadronDiffA");
    std::string outputFile = cl.Get("outputFile", "dNdEta_UPCDiffractivePythiaPrediction.pdf");

    TFile* file = TFile::Open(inputFile.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << inputFile << std::endl;
        return 1;
    }

    TH1D* hist = (TH1D*)(file->Get(histName));
    if (!hist) {
        std::cerr << "Error getting histogram: " << histName << std::endl;
        file->Close();
        return 1;
    }

    TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 600);
    TPad* pad = (TPad*) canvas->GetPad(0);
    pad->cd();
    SetTDRStyle();

    double marginLeft       = 0.16*1.2;
    double marginRight      = 0.04*1.2;
    double marginTop        = 0.06*1.2;
    double marginBottom     = 0.14*1.2;
    pad->SetMargin(marginLeft, marginRight, marginBottom, marginTop);

    hist->GetXaxis()->SetRangeUser(-14, 12);
    hist->GetXaxis()->SetTitle("#eta");
    SetTH1Fonts(hist, 0.8);
    hist->GetYaxis()->SetTitleOffset(1.6);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->SetLineColor(cmsBlue);
    hist->SetLineWidth(2);
    hist->GetXaxis()->SetNdivisions(5,0,5);

    hist->Draw("HIST");
    AddPythiaHeader_custom(pad, "", true, 0.06, -0.15);
    AddUPCHeader(pad, "9.62 TeV", "pO Diff. UPC", -1, 0.06);
    pad->Update();

    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(outputFile.c_str());

    file->Close();
    delete canvas;

    return 0;
}
