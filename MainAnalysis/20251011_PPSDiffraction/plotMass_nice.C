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

int main(int argc, char** argv) {
    CommandLine cl(argc, argv);

    std::string inputFile = cl.Get("inputFile", "output_loop_ppsTag_0_zdcPlusVeto_0_HFPlusGapVeto_1.root");
    TString histName = cl.Get("histName", "hMass");
    std::string outputFile = cl.Get("outputFile", "mass_plot.pdf");

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

    hist->GetXaxis()->SetRangeUser(0.97, 1.1);
    hist->GetXaxis()->SetTitle("#phi Mass (GeV)");
    SetTH1Fonts(hist, 0.75);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->SetLineColor(cmsBlue);
    hist->SetLineWidth(2);
    hist->GetXaxis()->SetNdivisions(5,0,5);

    hist->Draw("HIST");
    AddCMSHeader(pad, "Internal", true, 0.065);
    AddUPCHeader(pad, "9.62 TeV", "pO UPC", -1, 0.065);
    AddPlotLabel(pad, "|#eta_{#phi}| < 1.5", 0.065);
    pad->Update();

    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(outputFile.c_str());

    file->Close();
    delete canvas;

    return 0;
}