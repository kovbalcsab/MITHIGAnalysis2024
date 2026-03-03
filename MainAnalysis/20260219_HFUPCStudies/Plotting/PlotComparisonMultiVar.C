#include <TTree.h>
#include <TStyle.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>

#include "CommandLine.h" // Yi's Commandline bundle

void plotHistograms(std::vector<TH1D*> &hists, std::vector<Color_t> &colors, std::vector<std::string> &labels, std::string xTitle, std::string yTitle, std::string plotTitle, std::string outputFileName, 
    double MinDzeroPT, double MaxDzeroPT, double MinDzeroY, double MaxDzeroY, int IsGammaN,
    bool logY=false, double yMin=0, double yMax=0, double xMin=0, double xMax=0, bool normalize=false, int rebin=1, bool doRatio=false) {
    TCanvas* c1 = new TCanvas(Form("%s", hists[0]->GetName()), Form("%s", plotTitle.c_str()), 800, 800);

    if (hists.size() < 2) { doRatio = false; }

    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.04);
    c1->SetBottomMargin(0.12);
    c1->SetTopMargin(0.08);
    c1->cd();
    gStyle->SetOptStat(0);
    if (logY) { c1->SetLogy(); if (yMin <= 0) yMin = 1e-4; }

    TPad* padUp = nullptr;
    if (doRatio) {
        padUp = new TPad(Form("padUp%s", hists[0]->GetName()), "padUp", 0, 0.3, 1, 1);
        padUp->SetBottomMargin(0);
        padUp->Draw();
        padUp->cd();
        if (logY) padUp->SetLogy();
    }
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    auto styleUpHist = [](TH1 *hist) {
        hist->GetXaxis()->SetTitleSize(0.045);
        hist->GetXaxis()->SetLabelSize(0.04);
        hist->GetYaxis()->SetTitleSize(0.045);
        hist->GetYaxis()->SetLabelSize(0.04);
        hist->GetYaxis()->ChangeLabel(1,0,0);
    };

    auto styleDownHist = [](TH1 *hist) {
        hist->GetXaxis()->SetTitleSize(0.05 / 0.4);
        hist->GetXaxis()->SetLabelSize(0.045 / 0.4);
        hist->GetYaxis()->SetTitleSize(0.05 / 0.4);
        hist->GetYaxis()->SetLabelSize(0.045 / 0.4);
        hist->GetXaxis()->SetLabelOffset(0.01);
        hist->GetYaxis()->SetTitleOffset(.3);
        hist->GetYaxis()->ChangeLabel(-1,0,0);
        hist->GetYaxis()->SetNdivisions(505, true);
    };

    
    TLegend* legend1 = new TLegend(0.5, 0.65, 0.8, 0.85);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->SetTextSize(0.04);
    legend1->SetTextFont(42);

    for (size_t i = 0; i < hists.size(); ++i) {
        TH1D* currHist = hists[i];
        if (rebin > 1) {
            currHist->Rebin(rebin);
        }
        if (normalize && currHist->Integral("width") > 0) {
            currHist->Scale(1.0 / currHist->Integral("width")); // takes into account bin widths
        }
        currHist->SetLineColor(colors[i]);
        
        currHist->SetLineWidth(2);
        currHist->GetXaxis()->SetTitle(xTitle.c_str());
        currHist->GetXaxis()->SetRangeUser(xMin, xMax);
        currHist->GetYaxis()->SetRangeUser(yMin, yMax);
        currHist->GetYaxis()->SetTitle(yTitle.c_str());
        currHist->GetYaxis()->SetTitleOffset(1);
        currHist->SetStats(0);
        
        if (doRatio) {
            styleUpHist(currHist);
        }

        if (i == 0) {
            currHist->Draw("ep");
        } else {
            currHist->Draw("SAME ep");
        }
        legend1->AddEntry(currHist, labels[i].c_str(), "l");
    }
    legend1->Draw();
    // Make label for pt-y bin
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.SetNDC();
    if (MinDzeroPT > -999 && MaxDzeroPT > -999) latex.DrawLatex(0.15, 0.85, Form("%0.0f < D^{0} p_{T} < %0.0f GeV/c", MinDzeroPT, MaxDzeroPT));
    if (MinDzeroY > -999 && MaxDzeroY > -999) latex.DrawLatex(0.15, 0.80, Form("%0.2f < D^{0} y < %0.2f", MinDzeroY, MaxDzeroY));
    if (IsGammaN > -1) latex.DrawLatex(0.15, 0.75, Form("%s", IsGammaN ? "Gamma-N" : "N-Gamma"));

    TPad* padDown;
    if (doRatio) {
        c1->cd();
        padDown = new TPad(Form("padDown%s", hists[0]->GetName()), "padDown", 0, 0.05, 1, 0.3);
        padDown->SetTopMargin(0);
        padDown->SetBottomMargin(0.4);
        padDown->Draw();
        padDown->cd();
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);

        TH1D* denHist = (TH1D*)hists[0]->Clone(Form("ratio_%s", hists[0]->GetName()));
        for (size_t i = 1; i < hists.size(); ++i) {
            TH1D* numHist = hists[i];
            TH1D* ratioHist = (TH1D*)numHist->Clone(Form("ratio_%s_vs_%s", numHist->GetName(), denHist->GetName()));
            ratioHist->Divide(denHist);
            ratioHist->SetStats(0);
            styleDownHist(ratioHist);
            ratioHist->GetYaxis()->SetTitle("Ratio");
            ratioHist->GetYaxis()->SetRangeUser(0.0, 2);
            if (i == 1) {
                ratioHist->Draw("ep");
            } else {
                ratioHist->Draw("SAME ep");
            }
        }   

        // Draw unity line
        TH1D* unityLine = new TH1D("unityLine", "unityLine", 2, xMin-(xMax-xMin)/2, xMax+(xMax-xMin)/2);
        unityLine->SetBinContent(1, 1);
        unityLine->SetBinContent(2, 1);
        unityLine->SetLineColor(kGray+2);
        unityLine->SetLineStyle(2);
        unityLine->Draw("SAME L");
    }

    c1->Update();
    c1->SaveAs(outputFileName.c_str());
}

int main(int argc, char** argv) {
    CommandLine CL(argc, argv);

    std::vector<std::string> inputFileNames = CL.GetStringVector("InputFileNames", "");
    std::vector<std::string> varNames = CL.GetStringVector("VarNames", "");
    std::vector<bool> isTreeVar = CL.GetBoolVector("IsTreeVar", ""); // true if the variable is stored as a TTree, false if it's a TH1D
    std::vector<std::string> labels = CL.GetStringVector("Labels", "");
    // Replace ~ with space in labels
    for (auto& label : labels) {
        size_t pos = label.find('~');
        while (pos != std::string::npos) {
            label[pos] = ' ';
            pos = label.find('~');
        }
    }

    double MinDzeroPT = CL.GetDouble("MinDzeroPT", -999);
    double MaxDzeroPT = CL.GetDouble("MaxDzeroPT", -999);
    double MinDzeroY = CL.GetDouble("MinDzeroY", -999);
    double MaxDzeroY = CL.GetDouble("MaxDzeroY", -999);
    int IsGammaN = CL.GetInt("IsGammaN", -1); // -1 = not specified, 0 = N-Gamma, 1 = Gamma-N
    std::string outputFileName = CL.Get("OutputFileName", "CCF_output"); // without .root extension

    std::string xTitle = CL.Get("XTitle", "X-axis");
    std::string yTitle = CL.Get("YTitle", "Y-axis");
    std::string plotTitle = CL.Get("PlotTitle", "Comparison Plot");
    // Replace ~ with space in titles
    size_t pos = xTitle.find('~');
    while (pos != std::string::npos) {
        xTitle[pos] = ' ';
        pos = xTitle.find('~');
    }
    pos = yTitle.find('~');
    while (pos != std::string::npos) {
        yTitle[pos] = ' ';
        pos = yTitle.find('~');
    }
    pos = plotTitle.find('~');
    while (pos != std::string::npos) {
        plotTitle[pos] = ' ';
        pos = plotTitle.find('~');
    }

    bool logY = CL.GetBool("LogY", false);
    double yMin = CL.GetDouble("YMin", 0);
    double yMax = CL.GetDouble("YMax", 0);
    double xMin = CL.GetDouble("XMin", 0);
    double xMax = CL.GetDouble("XMax", 0);
    std::vector<int> nBins = CL.GetIntVector("NBins", ""); // Only neccessary if the input is TTree, otherwise the binning will be taken from the input histogram
    bool normalize = CL.GetBool("Normalize", false);
    int rebin = CL.GetInt("Rebin", 1);
    bool doRatio = CL.GetBool("DoRatio", true);
    
    std::vector<TH1D*> hists;
    std::vector<TFile*> inputFiles;
    std::vector<TTree*> inputTrees;
    std::vector<Color_t> colors = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange+7};

    // Load all input files and retrieve histograms
    for (int iFile = 0; iFile < inputFileNames.size(); ++iFile) {
        TFile* inFile = TFile::Open(inputFileNames[iFile].c_str(), "READ");
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error opening file: " << inputFileNames[iFile] << std::endl;
            continue;
        }

        if (isTreeVar[iFile]) {
            TTree* tree = dynamic_cast<TTree*>(inFile->Get("OutputTree"));
            if (!tree) {
                std::cerr << "Error: TTree " << varNames[iFile] << " not found in file " << inputFileNames[iFile] << std::endl;
                continue;
            }
            inputTrees.push_back(tree);
            // Create histogram from tree variable
            int nBinsForTree = (iFile < nBins.size()) ? nBins[iFile] : 400; // Default to 100 bins if not specified
            TH1D* histCurr = new TH1D(Form("hist_%d", iFile), "", nBinsForTree, xMin, xMax); // Adjust binning as needed
            histCurr->Sumw2();
            tree->Draw(Form("%s>>hist_%d", varNames[iFile].c_str(), iFile), "", "goff");
            hists.push_back(histCurr);
        } else {
            TH1D* histCurr = dynamic_cast<TH1D*>(inFile->Get(varNames[iFile].c_str()));
            if (histCurr) hists.push_back((TH1D*)histCurr->Clone(Form("hist_%d", iFile)));
            else {
                std::cerr << "Error: Histogram " << varNames[iFile] << " not found in file " << inputFileNames[iFile] << std::endl;
                continue;
            }
        }

        inputFiles.push_back(inFile);
    }

    if (hists.empty()) {
        std::cerr << "No histograms found in input files." << std::endl;
        return -1;
    }

    // Plot histograms
    plotHistograms(hists, colors, labels, xTitle, yTitle, plotTitle, outputFileName, 
        MinDzeroPT, MaxDzeroPT, MinDzeroY, MaxDzeroY, IsGammaN,
        logY, yMin, yMax, xMin, xMax, normalize, rebin, doRatio);

    for (auto* tree : inputTrees) {
        delete tree;
    }

    // Close input files
    for (auto* file : inputFiles) {
        file->Close();
    }

    return 0;
}
