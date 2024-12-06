#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLatex.h>

#include <iostream>
#include <string>

using namespace std;
#include "CommandLine.h" // Yi's Commandline bundle
#include "Messenger.h"   // Yi's Messengers for reading data files
#include "ProgressBar.h" // Yi's fish progress bar
#include "helpMessage.h" // Print out help message
#include "parameter.h"   // The parameters used in the analysis
#include "utilities.h"   // Yen-Jie's random utility functions


// function to plot prompt, non-prompt Dzero results
void plotDzero(string inputInclusive, string inputPrompt, string inputNonPrompt, string tag){
    // input files
    TFile* inclusiveFile = TFile::Open(inputInclusive.c_str());
    TFile* promptFile = TFile::Open(inputPrompt.c_str());
    TFile* nonpromptFile = TFile::Open(inputNonPrompt.c_str());

    // get the histogram
    TH1D* hInclusivePt = (TH1D*)inclusiveFile->Get("hDpt");
    hInclusivePt->SetName("hInclusivePt");
    TH1D* hInclusiveGenPt = (TH1D*)inclusiveFile->Get("hDptGen");
    hInclusiveGenPt->SetName("hInclusiveGenPt");
    TH1D* hPromptPt = (TH1D*)promptFile->Get("hDpt");
    hPromptPt->SetName("hPromptPt");
    TH1D* hPromptGenPt = (TH1D*)promptFile->Get("hDptGen");
    hPromptGenPt->SetName("hPromptGenPt");
    TH1D* hNonPromptPt = (TH1D*)nonpromptFile->Get("hDpt");
    hNonPromptPt->SetName("hNonPromptPt");
    TH1D* hNonPromptGenPt = (TH1D*)nonpromptFile->Get("hDptGen");
    hNonPromptGenPt->SetName("hNonPromptGenPt");

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);




    // ---------- reco histograms --------------------------------------------------
    TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->AddEntry(hInclusivePt, Form("Inclusive, Total Entries: %0.0f", hInclusivePt->GetEntries()));
    leg->AddEntry(hPromptPt, Form("Prompt, Total Entries: %0.0f", hPromptPt->GetEntries()));
    leg->AddEntry(hNonPromptPt, Form("Non-Prompt, Total Entries: %0.0f", hNonPromptPt->GetEntries()));

    // draw the histograms
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    c->SetLogy();
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetRightMargin(0.05);

    hInclusivePt->SetMarkerStyle(20);
    hInclusivePt->SetMarkerColor(kBlack);
    hInclusivePt->SetLineColor(kBlack);
    hPromptPt->SetMarkerStyle(20);
    hPromptPt->SetLineColor(kRed);
    hPromptPt->SetMarkerColor(kRed);
    hNonPromptPt->SetMarkerColor(kBlue);
    hNonPromptPt->SetLineColor(kBlue);
    hNonPromptPt->SetMarkerStyle(20);

    hInclusivePt->GetXaxis()->SetTitle("#it{D}^{0} #it{p}_{T} (GeV)");

    TLatex* cms = new TLatex(0.10,0.92,"#bf{CMS} #it{Simulation} PYTHIA 8 #gamma N Photon Beam B");
    cms->SetNDC();
    cms->SetTextSize(0.05);
    cms->SetTextFont(42);

    hInclusivePt->Draw();
    hPromptPt->Draw("same");
    hNonPromptPt->Draw("same");
    leg->Draw();
    cms->Draw();

    c->SaveAs(Form("pthat0_promptnonprompt_beamB_%s.pdf", tag.c_str()));

    // ---------- gen histograms --------------------------------------------------
    TLegend* legGen = new TLegend(0.5, 0.7, 0.9, 0.85);
    legGen->SetBorderSize(0);
    legGen->SetFillStyle(0);
    legGen->SetTextFont(42);
    legGen->AddEntry(hInclusivePt, Form("Inclusive, Total Entries: %0.0f", hInclusiveGenPt->GetEntries()));
    legGen->AddEntry(hPromptPt, Form("Prompt, Total Entries: %0.0f", hPromptGenPt->GetEntries()));
    legGen->AddEntry(hNonPromptPt, Form("Non-Prompt, Total Entries: %0.0f", hNonPromptGenPt->GetEntries()));

    // draw the histograms
    TCanvas* cGen = new TCanvas("cGen", "cGen", 800, 600);
    cGen->SetLogy();
    cGen->SetTickx(1);
    cGen->SetTicky(1);
    cGen->SetRightMargin(0.05);

    hInclusiveGenPt->SetMarkerStyle(20);
    hInclusiveGenPt->SetMarkerColor(kBlack);
    hInclusiveGenPt->SetLineColor(kBlack);
    hPromptGenPt->SetMarkerStyle(20);
    hPromptGenPt->SetLineColor(kRed);
    hPromptGenPt->SetMarkerColor(kRed);
    hNonPromptGenPt->SetMarkerColor(kBlue);
    hNonPromptGenPt->SetLineColor(kBlue);
    hNonPromptGenPt->SetMarkerStyle(20);

    hInclusiveGenPt->GetXaxis()->SetTitle("#it{D}^{0} Gen #it{p}_{T} (GeV)");

    hInclusiveGenPt->Draw();
    hPromptGenPt->Draw("same");
    hNonPromptGenPt->Draw("same");
    legGen->Draw();
    cms->Draw();

    cGen->SaveAs(Form("pthat0_promptnonprompt_beamB_gen_%s.pdf", tag.c_str()));


}

//============================================================//
// Main analysis
//============================================================//
int main(int argc, char *argv[]) {
  if (printHelpMessage(argc, argv))
    return 0;
  CommandLine CL(argc, argv);
  // input to this step is the result of running
  string inputPrompt = CL.Get("InputPrompt", "output_Prompt.root"); // Input file prompt
  string inputNonPrompt = CL.Get("InputNonPrompt", "output_NonPrompt.root"); // input file non-prompt
  string inputInclusive = CL.Get("Input", "output.root"); // inclusive input
  string tag = CL.Get("tag", "test");
  plotDzero(inputInclusive, inputPrompt, inputNonPrompt, tag);
}