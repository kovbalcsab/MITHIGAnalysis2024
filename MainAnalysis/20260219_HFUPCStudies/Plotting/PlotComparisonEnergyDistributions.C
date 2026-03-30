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
#include <TLine.h>
#include <TGraphErrors.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "CommandLine.h" // Yi's Commandline bundle

struct MetricMatrices
{
    TH2D *chi2PValue = nullptr;
    TH2D *wasserstein = nullptr;
    TH2D *jsDistance = nullptr;
    TH2D *hellinger = nullptr;
    TH2D *poissonDeviance = nullptr;
};

static bool compatibleBinning(const TH1D *a, const TH1D *b, double tolerance = 1e-12)
{
    if(a == nullptr || b == nullptr)
        return false;
    if(a->GetNbinsX() != b->GetNbinsX())
        return false;
    for(int edge = 1; edge <= a->GetNbinsX() + 1; ++edge)
    {
        const double aEdge = a->GetXaxis()->GetBinLowEdge(edge);
        const double bEdge = b->GetXaxis()->GetBinLowEdge(edge);
        if(std::fabs(aEdge - bEdge) > tolerance)
            return false;
    }
    return true;
}

static TH2D *makeLabeledMatrix(const std::string &name,
                               const std::string &title,
                               const std::vector<std::string> &labels)
{
    const int n = static_cast<int>(labels.size());
    TH2D *matrix = new TH2D(name.c_str(), title.c_str(), n, 0.5, n + 0.5, n, 0.5, n + 0.5);
    matrix->SetStats(0);
    for(int i = 0; i < n; ++i)
    {
        matrix->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
        matrix->GetYaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }
    return matrix;
}

static std::vector<double> buildProbabilityVector(const TH1D *hist,
                                                  double xMin, double xMax,
                                                  std::vector<double> *binWidths = nullptr,
                                                  double epsilon = 0.0)
{
    std::vector<double> masses;
    if(hist == nullptr)
        return masses;

    const bool restrictRange = (xMin < xMax);
    for(int b = 1; b <= hist->GetNbinsX(); ++b)
    {
        const double x = hist->GetBinCenter(b);
        if(restrictRange && (x < xMin || x > xMax))
            continue;
        const double width = hist->GetBinWidth(b);
        const double content = std::max(0.0, hist->GetBinContent(b));
        masses.push_back(content * width);
        if(binWidths != nullptr)
            binWidths->push_back(width);
    }
    if(masses.empty())
        return masses;

    double sum = 0.0;
    for(double m : masses)
        sum += m;
    if(sum <= 0.0)
        return std::vector<double>();

    for(double &m : masses)
        m /= sum;

    if(epsilon > 0.0)
    {
        double renorm = 0.0;
        for(double &m : masses)
        {
            m = std::max(m, epsilon);
            renorm += m;
        }
        if(renorm > 0.0)
        {
            for(double &m : masses)
                m /= renorm;
        }
    }

    return masses;
}

static double computeWasserstein1(const TH1D *a, const TH1D *b, double xMin, double xMax, double epsilon)
{
    if(!compatibleBinning(a, b))
        return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> widths;
    std::vector<double> p = buildProbabilityVector(a, xMin, xMax, &widths, epsilon);
    std::vector<double> q = buildProbabilityVector(b, xMin, xMax, nullptr, epsilon);
    if(p.empty() || q.empty() || p.size() != q.size() || p.size() != widths.size())
        return std::numeric_limits<double>::quiet_NaN();

    double cdfP = 0.0;
    double cdfQ = 0.0;
    double w1 = 0.0;
    for(size_t i = 0; i < p.size(); ++i)
    {
        cdfP += p[i];
        cdfQ += q[i];
        w1 += std::fabs(cdfP - cdfQ) * widths[i];
    }
    return w1;
}

static double computeJSDistance(const std::vector<double> &p, const std::vector<double> &q)
{
    if(p.size() != q.size() || p.empty())
        return std::numeric_limits<double>::quiet_NaN();

    double jsd = 0.0;
    for(size_t i = 0; i < p.size(); ++i)
    {
        const double pk = p[i];
        const double qk = q[i];
        const double mk = 0.5 * (pk + qk);
        if(pk > 0.0)
            jsd += 0.5 * pk * std::log(pk / mk);
        if(qk > 0.0)
            jsd += 0.5 * qk * std::log(qk / mk);
    }
    if(jsd < 0.0)
        jsd = 0.0;
    return std::sqrt(jsd);
}

static double computeHellinger(const std::vector<double> &p, const std::vector<double> &q)
{
    if(p.size() != q.size() || p.empty())
        return std::numeric_limits<double>::quiet_NaN();

    double sum = 0.0;
    for(size_t i = 0; i < p.size(); ++i)
    {
        const double diff = std::sqrt(std::max(0.0, p[i])) - std::sqrt(std::max(0.0, q[i]));
        sum += diff * diff;
    }
    return std::sqrt(0.5 * sum);
}

static double computePoissonDevianceOneWay(const TH1D *obs, const TH1D *ref,
                                           double xMin, double xMax,
                                           double floorMu)
{
    if(!compatibleBinning(obs, ref))
        return std::numeric_limits<double>::quiet_NaN();

    const bool restrictRange = (xMin < xMax);
    double dev = 0.0;
    for(int b = 1; b <= obs->GetNbinsX(); ++b)
    {
        const double x = obs->GetBinCenter(b);
        if(restrictRange && (x < xMin || x > xMax))
            continue;

        const double n = std::max(0.0, obs->GetBinContent(b));
        const double mu = std::max(floorMu, ref->GetBinContent(b));
        if(n > 0.0)
            dev += 2.0 * (mu - n + n * std::log(n / mu));
        else
            dev += 2.0 * mu;
    }
    return dev;
}

static double computeSymmetricPoissonDeviance(const TH1D *a, const TH1D *b,
                                              double xMin, double xMax,
                                              double floorMu)
{
    const double dAB = computePoissonDevianceOneWay(a, b, xMin, xMax, floorMu);
    const double dBA = computePoissonDevianceOneWay(b, a, xMin, xMax, floorMu);
    if(!std::isfinite(dAB) || !std::isfinite(dBA))
        return std::numeric_limits<double>::quiet_NaN();
    return 0.5 * (dAB + dBA);
}

static MetricMatrices buildPairwiseRobustMetricMatrices(const std::vector<TH1D*> &statAbsHists,
                                                        const std::vector<TH1D*> &statNormHists,
                                                        const std::vector<std::string> &labels,
                                                        double xMin, double xMax,
                                                        double epsilon,
                                                        double poissonFloorMu)
{
    MetricMatrices result;
    result.chi2PValue = makeLabeledMatrix("hChi2PValueMatrix",
                                          "Pairwise #chi^{2} p-value matrix;Distribution i;Distribution j", labels);
    result.wasserstein = makeLabeledMatrix("hWassersteinMatrix",
                                           "Pairwise Wasserstein-1 distance;Distribution i;Distribution j", labels);
    result.jsDistance = makeLabeledMatrix("hJSDistanceMatrix",
                                          "Pairwise Jensen-Shannon distance;Distribution i;Distribution j", labels);
    result.hellinger = makeLabeledMatrix("hHellingerMatrix",
                                         "Pairwise Hellinger distance;Distribution i;Distribution j", labels);
    result.poissonDeviance = makeLabeledMatrix("hPoissonDevianceMatrix",
                                               "Pairwise symmetric Poisson deviance;Distribution i;Distribution j", labels);

    const int n = static_cast<int>(labels.size());
    for(int i = 0; i < n; ++i)
    {
        for(int j = i; j < n; ++j)
        {
            double chi2P = 0.0;
            double w1 = 0.0;
            double js = 0.0;
            double h = 0.0;
            double pd = 0.0;

            if(i == j)
            {
                chi2P = 1.0;
                w1 = 0.0;
                js = 0.0;
                h = 0.0;
                pd = 0.0;
            }
            else
            {
                TH1D *absI = statAbsHists[i];
                TH1D *absJ = statAbsHists[j];
                TH1D *normI = statNormHists[i];
                TH1D *normJ = statNormHists[j];

                if(absI != nullptr && absJ != nullptr && normI != nullptr && normJ != nullptr
                   && compatibleBinning(absI, absJ) && compatibleBinning(normI, normJ))
                {
                    chi2P = normI->Chi2Test(normJ, "WW P");
                    w1 = computeWasserstein1(normI, normJ, xMin, xMax, epsilon);

                    std::vector<double> p = buildProbabilityVector(normI, xMin, xMax, nullptr, epsilon);
                    std::vector<double> q = buildProbabilityVector(normJ, xMin, xMax, nullptr, epsilon);
                    js = computeJSDistance(p, q);
                    h = computeHellinger(p, q);
                    pd = computeSymmetricPoissonDeviance(absI, absJ, xMin, xMax, poissonFloorMu);
                }
            }

            if(!std::isfinite(chi2P)) chi2P = 0.0;
            if(!std::isfinite(w1)) w1 = 0.0;
            if(!std::isfinite(js)) js = 0.0;
            if(!std::isfinite(h)) h = 0.0;
            if(!std::isfinite(pd)) pd = 0.0;

            result.chi2PValue->SetBinContent(i + 1, j + 1, chi2P);
            result.chi2PValue->SetBinContent(j + 1, i + 1, chi2P);

            result.wasserstein->SetBinContent(i + 1, j + 1, w1);
            result.wasserstein->SetBinContent(j + 1, i + 1, w1);

            result.jsDistance->SetBinContent(i + 1, j + 1, js);
            result.jsDistance->SetBinContent(j + 1, i + 1, js);

            result.hellinger->SetBinContent(i + 1, j + 1, h);
            result.hellinger->SetBinContent(j + 1, i + 1, h);

            result.poissonDeviance->SetBinContent(i + 1, j + 1, pd);
            result.poissonDeviance->SetBinContent(j + 1, i + 1, pd);
        }
    }

    return result;
}

static void plotPairwiseMatrixToPDF(TH2D *matrix, const std::string &outputFileName,
                                    double zMin, double zMax, bool autoMax = false)
{
    if(matrix == nullptr)
        return;

    double minUse = zMin;
    double maxUse = zMax;
    if(autoMax)
    {
        maxUse = 0.0;
        for(int ix = 1; ix <= matrix->GetNbinsX(); ++ix)
            for(int iy = 1; iy <= matrix->GetNbinsY(); ++iy)
                maxUse = std::max(maxUse, matrix->GetBinContent(ix, iy));
        if(maxUse <= minUse)
            maxUse = minUse + 1.0;
    }

    TCanvas *c = new TCanvas(Form("cStatMatrix_%s", matrix->GetName()), "Pairwise statistical matrix", 900, 800);
    c->SetLeftMargin(0.16);
    c->SetRightMargin(0.16);
    c->SetBottomMargin(0.16);
    c->SetTopMargin(0.08);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".3f");
    matrix->SetMinimum(minUse);
    matrix->SetMaximum(maxUse);
    matrix->Draw("COLZ TEXT");
    c->SaveAs(outputFileName.c_str());
    delete c;
}

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
    if (logY) {
        if (!doRatio) c1->SetLogy();
        if (yMin <= 0) yMin = 1e-4;
    }

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
        const Color_t color = colors.empty() ? kBlack : colors[i % colors.size()];
        currHist->SetLineColor(color);
        currHist->SetMarkerColor(color);
        currHist->SetMarkerStyle(20 + (i % 5));
        currHist->SetMarkerSize(0.95);
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
            currHist->Draw("HIST");
            currHist->Draw("SAME E1P");
        } else {
            currHist->Draw("SAME HIST");
            currHist->Draw("SAME E1P");
        }
        legend1->AddEntry(currHist, labels[i].c_str(), "lp");
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

    TPad* padDown = nullptr;
    TH1D *denHistToDelete = nullptr;
    TH1D *axisHistToDelete = nullptr;
    TLine *lineUnityToDelete = nullptr;
    std::vector<TGraphErrors*> ratioGraphsToDelete;
    if (doRatio) {
        c1->cd();
        padDown = new TPad(Form("padDown%s", hists[0]->GetName()), "padDown", 0, 0.05, 1, 0.3);
        padDown->SetTopMargin(0);
        padDown->SetBottomMargin(0.4);
        padDown->SetLogy(false);
        padDown->Draw();
        padDown->cd();
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);

        TH1D* denHist = (TH1D*)hists[0]->Clone(Form("ratio_den_%s", hists[0]->GetName()));
        denHist->SetDirectory(nullptr);
        denHistToDelete = denHist;
        ratioGraphsToDelete.reserve(hists.size());

        TH1D* axisHist = (TH1D*)hists[0]->Clone(Form("ratio_axis_%s", hists[0]->GetName()));
        axisHist->Reset("ICES");
        axisHist->SetDirectory(nullptr);
        axisHistToDelete = axisHist;
        axisHist->SetStats(0);
        styleDownHist(axisHist);
        axisHist->GetYaxis()->SetTitle("Ratio to first");
        axisHist->GetYaxis()->SetRangeUser(0.0, 2.0);
        axisHist->GetXaxis()->SetTitle(xTitle.c_str());
        axisHist->GetXaxis()->SetRangeUser(xMin, xMax);
        axisHist->Draw("AXIS");

        TLine *lineUnity = new TLine(xMin, 1.0, xMax, 1.0);
        lineUnityToDelete = lineUnity;
        lineUnity->SetLineColor(kGray + 2);
        lineUnity->SetLineStyle(2);
        lineUnity->Draw("SAME");

        for (size_t i = 1; i < hists.size(); ++i) {
            TH1D* numHist = hists[i];
            TGraphErrors *gr = new TGraphErrors();
            gr->SetName(Form("gRatio_%zu", i));
            const Color_t color = colors.empty() ? kBlack : colors[i % colors.size()];
            gr->SetLineColor(color);
            gr->SetMarkerColor(color);
            gr->SetMarkerStyle(20 + (i % 5));
            gr->SetLineWidth(2);

            int p = 0;
            for(int b = 1; b <= numHist->GetNbinsX(); ++b)
            {
                const double x = numHist->GetBinCenter(b);
                if(x < xMin || x > xMax)
                    continue;
                const double ex = numHist->GetBinWidth(b) / 2.0;
                const double num = numHist->GetBinContent(b);
                const double en = numHist->GetBinError(b);
                int denBin = denHist->GetXaxis()->FindFixBin(x);
                if(denBin < 1)
                    denBin = 1;
                if(denBin > denHist->GetNbinsX())
                    denBin = denHist->GetNbinsX();
                const double den = denHist->GetBinContent(denBin);
                const double ed = denHist->GetBinError(denBin);
                if(den <= 0)
                    continue;
                const double ratio = num / den;
                double eratio = 0.0;
                if(num > 0)
                    eratio = ratio * std::sqrt((en/num)*(en/num) + (ed/den)*(ed/den));
                else
                    eratio = ratio * (ed/den);
                gr->SetPoint(p, x, ratio);
                gr->SetPointError(p, ex, eratio);
                p++;
            }
            gr->Draw("EP SAME");
            ratioGraphsToDelete.push_back(gr);
        }
    }

    c1->Update();
    c1->SaveAs(outputFileName.c_str());

    for(auto *g : ratioGraphsToDelete) delete g;
    delete lineUnityToDelete;
    delete axisHistToDelete;
    delete denHistToDelete;
    delete padDown;
    delete padUp;
    delete c1;
}

int main(int argc, char** argv) {
    CommandLine CL(argc, argv);

    std::vector<std::string> inputFileNames = CL.GetStringVector("InputFileNames", "");
    std::vector<std::string> varNames = CL.GetStringVector("VarNames", "");
    std::vector<bool> isTreeVar = CL.GetBoolVector("IsTreeVar", ""); // true if the variable is stored as a TTree, false if it's a TH1D
    std::vector<std::string> labels = CL.GetStringVector("Labels", "");
    if(inputFileNames.size() != varNames.size() || inputFileNames.size() != isTreeVar.size())
    {
        std::cerr << "InputFileNames, VarNames and IsTreeVar must have the same size." << std::endl;
        return -1;
    }
    if(!labels.empty() && labels.size() != inputFileNames.size())
    {
        std::cerr << "Labels must be empty or have the same size as InputFileNames." << std::endl;
        return -1;
    }
    if(labels.empty())
    {
        labels.resize(inputFileNames.size());
        for(size_t i = 0; i < labels.size(); ++i)
            labels[i] = Form("Dist %zu", i + 1);
    }
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
    std::string statMatrixOutputFileName = CL.Get("StatMatrixOutputFileName", outputFileName + "_Chi2PValueMatrix.pdf");
    std::string wassersteinOutputFileName = CL.Get("WassersteinOutputFileName", outputFileName + "_WassersteinMatrix.pdf");
    std::string jsDistanceOutputFileName = CL.Get("JSDistanceOutputFileName", outputFileName + "_JSDistanceMatrix.pdf");
    std::string hellingerOutputFileName = CL.Get("HellingerOutputFileName", outputFileName + "_HellingerMatrix.pdf");
    std::string poissonDevianceOutputFileName = CL.Get("PoissonDevianceOutputFileName", outputFileName + "_PoissonDevianceMatrix.pdf");
    int statRebin = CL.GetInt("StatRebin", 1);
    if(statRebin < 1) statRebin = 1;
    double statEpsilon = CL.GetDouble("StatEpsilon", 1e-12);
    if(statEpsilon < 0.0) statEpsilon = 0.0;
    double poissonFloorMu = CL.GetDouble("PoissonFloorMu", 1e-9);
    if(poissonFloorMu <= 0.0) poissonFloorMu = 1e-9;
    
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
            int nBinsForTree = (iFile < nBins.size()) ? nBins[iFile] : 100; // Default to 100 bins if not specified
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

    // Build stat clones before plotting mutations
    std::vector<TH1D*> statAbsHists;
    std::vector<TH1D*> statNormHists;
    statAbsHists.reserve(hists.size());
    statNormHists.reserve(hists.size());
    for(size_t i = 0; i < hists.size(); ++i)
    {
        TH1D *absClone = dynamic_cast<TH1D*>(hists[i]->Clone(Form("hStatAbs_%zu", i)));
        TH1D *normClone = dynamic_cast<TH1D*>(hists[i]->Clone(Form("hStatNorm_%zu", i)));
        if(absClone == nullptr || normClone == nullptr)
        {
            std::cerr << "Failed to clone histogram for statistic computation at index " << i << std::endl;
            delete absClone;
            delete normClone;
            for(auto *h : statAbsHists) delete h;
            for(auto *h : statNormHists) delete h;
            for (auto* file : inputFiles) file->Close();
            return -1;
        }
        absClone->SetDirectory(nullptr);
        normClone->SetDirectory(nullptr);
        if(statRebin > 1)
        {
            absClone->Rebin(statRebin);
            normClone->Rebin(statRebin);
        }
        const double integral = normClone->Integral("width");
        if(integral > 0.0)
            normClone->Scale(1.0 / integral);
        statAbsHists.push_back(absClone);
        statNormHists.push_back(normClone);
    }

    // Plot histograms
    plotHistograms(hists, colors, labels, xTitle, yTitle, plotTitle, outputFileName, 
        MinDzeroPT, MaxDzeroPT, MinDzeroY, MaxDzeroY, IsGammaN,
        logY, yMin, yMax, xMin, xMax, normalize, rebin, doRatio);

    MetricMatrices metricMatrices = buildPairwiseRobustMetricMatrices(
        statAbsHists, statNormHists, labels, xMin, xMax, statEpsilon, poissonFloorMu);

    plotPairwiseMatrixToPDF(metricMatrices.chi2PValue, statMatrixOutputFileName, 0.0, 1.0, false);
    plotPairwiseMatrixToPDF(metricMatrices.wasserstein, wassersteinOutputFileName, 0.0, 0.0, true);
    plotPairwiseMatrixToPDF(metricMatrices.jsDistance, jsDistanceOutputFileName, 0.0, 1.0, false);
    plotPairwiseMatrixToPDF(metricMatrices.hellinger, hellingerOutputFileName, 0.0, 1.0, false);
    plotPairwiseMatrixToPDF(metricMatrices.poissonDeviance, poissonDevianceOutputFileName, 0.0, 0.0, true);

    delete metricMatrices.chi2PValue;
    delete metricMatrices.wasserstein;
    delete metricMatrices.jsDistance;
    delete metricMatrices.hellinger;
    delete metricMatrices.poissonDeviance;
    for(auto *h : statAbsHists) delete h;
    for(auto *h : statNormHists) delete h;

    // Close input files
    for (auto* file : inputFiles) {
        file->Close();
    }

    return 0;
}
