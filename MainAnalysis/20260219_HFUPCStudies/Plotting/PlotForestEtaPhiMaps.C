#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TMath.h"
#include "TF1.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "CommandLine.h"

using namespace std;

////////////////////////////////////////////////////////////
// CMS style (minimal clean version)
////////////////////////////////////////////////////////////
void SetCMSStyle()
{
    gStyle->SetOptStat(0);

    gStyle->SetTitleSize(0.05,"XYZ");
    gStyle->SetLabelSize(0.045,"XYZ");

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetFrameLineWidth(1);

    gStyle->SetEndErrorSize(0);
}

////////////////////////////////////////////////////////////
// Histogram styling
////////////////////////////////////////////////////////////
void StyleHistogram(TH1D* hist, int color)
{
    if(!hist) return;

    hist->SetLineColor(color);
    hist->SetLineWidth(1);

    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.8);

    hist->GetXaxis()->SetTitle("HF E_{max} (GeV)");
    hist->GetYaxis()->SetTitle("Entries");

    hist->SetTitle("");

    hist->GetXaxis()->SetRangeUser(0,20);

    hist->GetYaxis()->SetRangeUser(1e-5, 1.);
}

////////////////////////////////////////////////////////////
// Draw histogram with smart axis visibility
////////////////////////////////////////////////////////////
void DrawHistogram(TH1D* hist,
                   bool showX,
                   bool showY,
                   int color)
{
    StyleHistogram(hist,color);

    if(!showX)
    {
        hist->GetXaxis()->SetLabelSize(0);
        hist->GetXaxis()->SetTitleSize(0);
    }

    if(!showY)
    {
        hist->GetYaxis()->SetLabelSize(0);
        hist->GetYaxis()->SetTitleSize(0);
    }

    hist->Draw("hist");
}

////////////////////////////////////////////////////////////

double gamma_pdf(double x, double k, double theta) {
    if (x < 0) return 0;
    return pow(x, k-1) * exp(-x/theta) / (pow(theta, k) * TMath::Gamma(k));
}

Double_t exponential_pdf(Double_t *x, Double_t *par) {
    if (x[0] < 0) return 0;
    return par[1]/par[0] * exp(-x[0]/par[0]);
}

Double_t fitf(Double_t *x,Double_t *par){
    return par[3]*gamma_pdf(x[0]-par[5], par[0], par[1])+par[4]*exp(-x[0]/par[2])/par[2];
}


int main(int argc, char *argv[])
{
    CommandLine CL(argc,argv);

    string InputFileName  = CL.Get("Input");
    string OutputFileName = CL.Get("Output");
    bool doCumulative = CL.GetBool("DoCumulative",false);

    SetCMSStyle();

    TFile *InputFile = TFile::Open(InputFileName.c_str());

    TF1* fitFunc = new TF1("fitFunc", fitf, 1, 20, 6);
    fitFunc->SetParNames("k", "theta", "lambda", "A_gamma", "A_exp", "offset");
    fitFunc->SetParameters(1, .3, 2, 0.05, 0.005, 0.);
    fitFunc->SetParLimits(0, 0.01, 20);
    fitFunc->SetParLimits(1, 0.1, 10);
    fitFunc->SetParLimits(2, 2, 10);
    fitFunc->SetParLimits(3, 0, 10);
    fitFunc->SetParLimits(4, 0, 10);
    fitFunc->SetParLimits(5, -1, 3);

    TF1* fitFuncExp = new TF1("fitFuncExp", exponential_pdf, 8, 20, 2);
    fitFuncExp->SetParameters(2, 0.05);
    fitFuncExp->SetParNames("lambda", "A_exp");
    fitFuncExp->SetParLimits(0, 2, 10);
    fitFuncExp->SetParLimits(1, 0, 10);

    ////////////////////////////////////////////////////////////
    // binning
    ////////////////////////////////////////////////////////////

    vector<float> etaBorders =
    {3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2};

    vector<float> phiBorders =
    {-M_PI,-2*M_PI/3,-M_PI/3,0,M_PI/3,2*M_PI/3,M_PI};

    const int nEta = etaBorders.size()-1;
    const int nPhi = phiBorders.size()-1;

    ////////////////////////////////////////////////////////////
    // load histograms
    ////////////////////////////////////////////////////////////

    vector<vector<TH1D*>> hPlus(nEta,vector<TH1D*>(nPhi,nullptr));
    vector<vector<TH1D*>> hMinus(nEta,vector<TH1D*>(nPhi,nullptr));
    vector<vector<TH1D*>> hPlusCum(nEta,vector<TH1D*>(nPhi,nullptr));
    vector<vector<TH1D*>> hMinusCum(nEta,vector<TH1D*>(nPhi,nullptr));

    for(int iEta=0;iEta<nEta;iEta++)
    for(int iPhi=0;iPhi<nPhi;iPhi++)
    {
        hPlus[iEta][iPhi] =
        (TH1D*)InputFile->Get(
            Form("hHFEMaxPlus_eta%.1f_%.1f_phi%.1f_%.1f",
            etaBorders[iEta],etaBorders[iEta+1],
            phiBorders[iPhi],phiBorders[iPhi+1]));  
            
        if(hPlus[iEta][iPhi]->Integral() > 0)
            hPlus[iEta][iPhi]->Scale(1.0 / hPlus[iEta][iPhi]->Integral());
        if (doCumulative)
            hPlusCum[iEta][iPhi]= (TH1D*)hPlus[iEta][iPhi]->GetCumulative(kFALSE);
    
        hMinus[iEta][iPhi] =
        (TH1D*)InputFile->Get(
            Form("hHFEMaxMinus_eta%.1f_%.1f_phi%.1f_%.1f",
            -etaBorders[iEta+1],-etaBorders[iEta],
            phiBorders[iPhi],phiBorders[iPhi+1]));
        if(hMinus[iEta][iPhi]->Integral() > 0)
            hMinus[iEta][iPhi]->Scale(1.0 / hMinus[iEta][iPhi]->Integral());
        if (doCumulative)
            hMinusCum[iEta][iPhi]= (TH1D*)hMinus[iEta][iPhi]->GetCumulative(kFALSE);
    }

    ////////////////////////////////////////////////////////////
    // pad geometry (fix squishing)
    ////////////////////////////////////////////////////////////

    int padW = 300;
    int padH = 300;

    double leftMargin   = 0.16;
    double rightMargin  = 0.02;
    double bottomMargin = 0.16;
    double topMargin    = 0.02;

    int canvasW = padW*nEta/(1-leftMargin-rightMargin);
    int canvasH = padH*nPhi/(1-bottomMargin-topMargin);

    TCanvas *cPlus  = new TCanvas("cPlus","",canvasW,canvasH);
    TCanvas *cMinus = new TCanvas("cMinus","",canvasW,canvasH);

    cPlus->Divide(nEta,nPhi,0,0);
    cMinus->Divide(nEta,nPhi,0,0);

    TLatex *label = new TLatex();
    label->SetNDC();
    label->SetTextSize(0.05);
    label->SetTextFont(42);

    std::vector<std::vector<TF1*>> fitResultsPlus(nEta, std::vector<TF1*>(nPhi, nullptr));
    std::vector<std::vector<TF1*>> fitResultsMinus(nEta, std::vector<TF1*>(nPhi, nullptr));

    ////////////////////////////////////////////////////////////
    // drawing
    ////////////////////////////////////////////////////////////

    for(int iEta=0;iEta<nEta;iEta++)
    for(int iPhi=0;iPhi<nPhi;iPhi++)
    {
        fitFunc->SetParameters(1, .3, 2, 0.05, 0.005, 0.);

        int padID = iPhi*nEta + iEta + 1;

        bool showY = (iEta==0);
        bool showX = (iPhi==nPhi-1);

        ///////////////////////
        // PLUS
        ///////////////////////

        cPlus->cd(padID);

        gPad->SetLogy();

        gPad->SetLeftMargin(showY ? leftMargin : 0.02);
        gPad->SetRightMargin(rightMargin);
        gPad->SetBottomMargin(showX ? bottomMargin : 0.02);
        gPad->SetTopMargin(topMargin);

        gPad->SetFrameLineWidth(1);

        TH1D* histToDrawPlus = doCumulative ? hPlusCum[iEta][iPhi] : hPlus[iEta][iPhi];
        DrawHistogram(histToDrawPlus,showX,showY,kBlue);
        //hPlus[iEta][iPhi]->Fit(fitFunc,"R");
        //fitResultsPlus[iEta][iPhi] = (TF1*)fitFunc->Clone();
        //fitResultsPlus[iEta][iPhi]->Draw("same");
        label->DrawLatex(0.15,0.85,
            Form("#eta: [%.1f,%.1f], #phi: [%.1f,%.1f]",
            etaBorders[iEta],etaBorders[iEta+1],
            phiBorders[iPhi],phiBorders[iPhi+1]));

        ///////////////////////
        // MINUS
        ///////////////////////

        fitFunc->SetParameters(1, .3, 2, 0.05, 0.005, 0.);
        cMinus->cd(padID);

        gPad->SetLogy();

        gPad->SetLeftMargin(showY ? leftMargin : 0.02);
        gPad->SetRightMargin(rightMargin);
        gPad->SetBottomMargin(showX ? bottomMargin : 0.02);
        gPad->SetTopMargin(topMargin);

        gPad->SetFrameLineWidth(1);

        TH1D* histToDrawMinus = doCumulative ? hMinusCum[iEta][iPhi] : hMinus[iEta][iPhi];
        DrawHistogram(histToDrawMinus,showX,showY,kRed);
        //hMinus[iEta][iPhi]->Fit(fitFunc,"R");
        //fitResultsMinus[iEta][iPhi] = (TF1*)fitFunc->Clone();
        //fitResultsMinus[iEta][iPhi]->Draw("same");
        label->DrawLatex(0.15,0.85,
            Form("#eta: [%.1f,%.1f], #phi: [%.1f,%.1f]",
            -etaBorders[iEta+1],-etaBorders[iEta],
            phiBorders[iPhi],phiBorders[iPhi+1]));
    }

    ////////////////////////////////////////////////////////////
    // save
    ////////////////////////////////////////////////////////////

    cPlus->SaveAs(
        Form("Plotting/ForestEtaPhiMapsPlots/HFEMaxPlusMaps_%s.pdf",
        OutputFileName.c_str()));

    cMinus->SaveAs(
        Form("Plotting/ForestEtaPhiMapsPlots/HFEMaxMinusMaps_%s.pdf",
        OutputFileName.c_str()));

    cPlus->SaveAs(
        Form("Plotting/ForestEtaPhiMapsPlots/HFEMaxPlusMaps_%s.png",
        OutputFileName.c_str()));

    cMinus->SaveAs(
        Form("Plotting/ForestEtaPhiMapsPlots/HFEMaxMinusMaps_%s.png",
        OutputFileName.c_str()));

    return 0;
}