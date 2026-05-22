#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "CommandLine.h"

namespace
{
   struct DistributionSet
   {
      TH2D *H36 = nullptr;
      TH2D *H18Main = nullptr;
      TH2D *H18WrapLow = nullptr;
      TH2D *H18WrapHigh = nullptr;
   };

   struct DisplayRange
   {
      double Minimum = 0.0;
      double Maximum = 1.0;
      bool UseLogZ = false;
   };

   bool PrintHelpMessage(int argc, char *argv[])
   {
      for(int i = 1; i < argc; ++i)
      {
         if(std::string(argv[i]) != "--help")
            continue;

         std::cout
            << "Usage: " << argv[0] << " [OPTIONS]\n"
            << "Options:\n"
            << "  --Input <path>               Input ROOT file. Default: output/HFDistributions_First_100.root\n"
            << "  --OutputDir <path>           Output directory. Default: output\n"
            << "  --OutputPrefix <prefix>      Output file prefix. Default: HFDistributions\n"
            << "  --OutputFormat <pdf|png>     Output image format. Default: pdf\n"
            << "  --HistogramSuffix <suffix>   String appended to expected histogram names. Default: empty\n"
            << "  --IgnoreEmptyWrapHigh <bool> If true, drop WrapHigh only when it is completely empty. Default: false\n"
            ;
         return true;
      }
      return false;
   }

   TH2D *GetHistogram(TFile &InputFile, const std::string &Name, bool PrintError = true)
   {
      TH2D *Histogram = dynamic_cast<TH2D *>(InputFile.Get(Name.c_str()));
      if(Histogram == nullptr)
      {
         if(PrintError == true)
            std::cerr << "Failed to load histogram " << Name << std::endl;
         return nullptr;
      }

      Histogram->SetDirectory(nullptr);
      return Histogram;
   }

   void DeleteSet(DistributionSet &Set)
   {
      delete Set.H36;
      delete Set.H18Main;
      delete Set.H18WrapLow;
      delete Set.H18WrapHigh;
      Set = {};
   }

   bool IsCompletelyEmpty(const TH2D *Histogram)
   {
      if(Histogram == nullptr)
         return true;

      for(int ix = 1; ix <= Histogram->GetNbinsX(); ++ix)
      {
         for(int iy = 1; iy <= Histogram->GetNbinsY(); ++iy)
         {
            if(Histogram->GetBinContent(ix, iy) != 0.0)
               return false;
         }
      }

      return true;
   }

   bool LoadSet(TFile &InputFile, const std::string &SidePrefix, const std::string &Suffix,
      const std::string &HistogramSuffix, DistributionSet &Set, bool IgnoreEmptyWrapHigh)
   {
      Set.H36 = GetHistogram(InputFile, "h36" + SidePrefix + Suffix + HistogramSuffix, false);
      Set.H18Main = GetHistogram(InputFile, "h18" + SidePrefix + "Main" + Suffix + HistogramSuffix, false);
      Set.H18WrapLow = GetHistogram(InputFile, "h18" + SidePrefix + "WrapLow" + Suffix + HistogramSuffix, false);
      Set.H18WrapHigh = GetHistogram(InputFile, "h18" + SidePrefix + "WrapHigh" + Suffix + HistogramSuffix, false);

      if(Set.H36 == nullptr || Set.H18Main == nullptr || Set.H18WrapLow == nullptr || Set.H18WrapHigh == nullptr)
      {
         DeleteSet(Set);
         return false;
      }

      if(IgnoreEmptyWrapHigh == true && IsCompletelyEmpty(Set.H18WrapHigh) == true)
      {
         delete Set.H18WrapHigh;
         Set.H18WrapHigh = nullptr;
      }

      return true;
   }

   void StyleHistogram(TH2D *Histogram, const std::string &Title)
   {
      Histogram->SetTitle(Title.c_str());
      Histogram->GetXaxis()->SetTitle("#eta");
      Histogram->GetYaxis()->SetTitle("#phi");
      Histogram->GetYaxis()->SetTitleOffset(1.1);
   }

   std::vector<TH2D *> AsVector(const DistributionSet &Set)
   {
      std::vector<TH2D *> Histograms;
      if(Set.H36 != nullptr)
         Histograms.push_back(Set.H36);
      if(Set.H18Main != nullptr)
         Histograms.push_back(Set.H18Main);
      if(Set.H18WrapLow != nullptr)
         Histograms.push_back(Set.H18WrapLow);
      if(Set.H18WrapHigh != nullptr)
         Histograms.push_back(Set.H18WrapHigh);
      return Histograms;
   }

   double GetAxisMinimum(const std::vector<TH2D *> &Histograms, bool UseX)
   {
      double Minimum = 0.0;
      bool Initialized = false;

      for(TH2D *Histogram : Histograms)
      {
         const TAxis *Axis = (UseX == true) ? Histogram->GetXaxis() : Histogram->GetYaxis();
         const double LowEdge = Axis->GetBinLowEdge(1);
         if(Initialized == false || LowEdge < Minimum)
            Minimum = LowEdge;
         Initialized = true;
      }

      return Minimum;
   }

   double GetAxisMaximum(const std::vector<TH2D *> &Histograms, bool UseX)
   {
      double Maximum = 0.0;
      bool Initialized = false;

      for(TH2D *Histogram : Histograms)
      {
         const TAxis *Axis = (UseX == true) ? Histogram->GetXaxis() : Histogram->GetYaxis();
         const double HighEdge = Axis->GetBinUpEdge(Axis->GetNbins());
         if(Initialized == false || HighEdge > Maximum)
            Maximum = HighEdge;
         Initialized = true;
      }

      return Maximum;
   }

   double GetHistogramMaximum(const std::vector<TH2D *> &Histograms)
   {
      double MaximumContent = 0.0;

      for(TH2D *Histogram : Histograms)
         MaximumContent = std::max(MaximumContent, Histogram->GetMaximum());

      return MaximumContent;
   }

   std::vector<double> GetNonZeroContents(const std::vector<TH2D *> &Histograms)
   {
      std::vector<double> Values;

      for(TH2D *Histogram : Histograms)
      {
         for(int ix = 1; ix <= Histogram->GetNbinsX(); ++ix)
         {
            for(int iy = 1; iy <= Histogram->GetNbinsY(); ++iy)
            {
               const double Value = Histogram->GetBinContent(ix, iy);
               if(Value > 0.0)
                  Values.push_back(Value);
            }
         }
      }

      std::sort(Values.begin(), Values.end());
      return Values;
   }

   double GetQuantile(const std::vector<double> &Values, double Fraction)
   {
      if(Values.empty() == true)
         return 0.0;

      const double Clamped = std::max(0.0, std::min(1.0, Fraction));
      const size_t Index = static_cast<size_t>(Clamped * static_cast<double>(Values.size() - 1));
      return Values[Index];
   }

   DisplayRange DetermineDisplayRange(const std::vector<TH2D *> &Histograms, const std::string &ModeLabel)
   {
      const std::vector<double> NonZeroValues = GetNonZeroContents(Histograms);
      const double Maximum = GetHistogramMaximum(Histograms);

      DisplayRange Result;
      if(NonZeroValues.empty() == true || Maximum <= 0.0)
      {
         Result.Minimum = 0.0;
         Result.Maximum = 1.0;
         Result.UseLogZ = false;
         return Result;
      }

      if(ModeLabel == "CountCrossCheck")
      {
         Result.Minimum = std::max(1.0, GetQuantile(NonZeroValues, 0.01));
         Result.Maximum = GetQuantile(NonZeroValues, 0.99);
         Result.UseLogZ = true;
      }
      else
      {
         Result.Minimum = GetQuantile(NonZeroValues, 0.01);
         Result.Maximum = GetQuantile(NonZeroValues, 0.99);
         Result.UseLogZ = false;
      }

      if(Result.Maximum <= Result.Minimum)
         Result.Maximum = Maximum;
      if(Result.Maximum <= Result.Minimum)
         Result.Maximum = Result.Minimum * 10.0;

      return Result;
   }

   TH2D BuildPaletteHelper(const std::vector<TH2D *> &Histograms,
      const std::string &Title,
      double XMinimum, double XMaximum,
      double YMinimum, double YMaximum,
      double FillValue)
   {
      TH2D Helper("PaletteFrame", Title.c_str(), 800, XMinimum, XMaximum, 800, YMinimum, YMaximum);
      Helper.SetDirectory(nullptr);
      StyleHistogram(&Helper, Title);
      Helper.SetStats(0);

      for(TH2D *Histogram : Histograms)
      {
         for(int ix = 1; ix <= Histogram->GetNbinsX(); ++ix)
         {
            const double XLow = Histogram->GetXaxis()->GetBinLowEdge(ix);
            const double XHigh = Histogram->GetXaxis()->GetBinUpEdge(ix);
            const int XBinMin = Helper.GetXaxis()->FindBin(XLow + 1e-9);
            const int XBinMax = Helper.GetXaxis()->FindBin(XHigh - 1e-9);

            for(int iy = 1; iy <= Histogram->GetNbinsY(); ++iy)
            {
               const double YLow = Histogram->GetYaxis()->GetBinLowEdge(iy);
               const double YHigh = Histogram->GetYaxis()->GetBinUpEdge(iy);
               const int YBinMin = Helper.GetYaxis()->FindBin(YLow + 1e-9);
               const int YBinMax = Helper.GetYaxis()->FindBin(YHigh - 1e-9);

               for(int dx = XBinMin; dx <= XBinMax; ++dx)
               {
                  for(int dy = YBinMin; dy <= YBinMax; ++dy)
                     Helper.SetBinContent(dx, dy, FillValue);
               }
            }
         }
      }

      return Helper;
   }

   void DrawSidePDF(const std::string &OutputPDF,
      const std::string &CanvasTitle,
      const DistributionSet &Set,
      const std::string &ModeLabel)
   {
      std::vector<TH2D *> Histograms = AsVector(Set);
      const double XMinimum = GetAxisMinimum(Histograms, true);
      const double XMaximum = GetAxisMaximum(Histograms, true);
      const double YMinimum = GetAxisMinimum(Histograms, false);
      const double YMaximum = GetAxisMaximum(Histograms, false);
      const DisplayRange Range = DetermineDisplayRange(Histograms, ModeLabel);

      TH2D PaletteFrame = BuildPaletteHelper(Histograms, CanvasTitle, XMinimum, XMaximum, YMinimum, YMaximum, Range.Minimum);
      PaletteFrame.SetMinimum(Range.Minimum);
      PaletteFrame.SetMaximum(Range.Maximum);

      TCanvas Canvas("Canvas", "", 1100, 900);
      Canvas.cd();
      gPad->SetRightMargin(0.16);
      gPad->SetLeftMargin(0.12);
      gPad->SetBottomMargin(0.12);
      gPad->SetLogz(Range.UseLogZ == true ? 1 : 0);
      PaletteFrame.Draw("COLZ");

      for(size_t i = 0; i < Histograms.size(); ++i)
      {
         Histograms[i]->SetMinimum(Range.Minimum);
         Histograms[i]->SetMaximum(Range.Maximum);

         Histograms[i]->Draw("COL SAME");
      }

      gPad->RedrawAxis();

      Canvas.Print(OutputPDF.c_str());
      std::cout << "Wrote " << CanvasTitle << " output: " << OutputPDF << std::endl;
   }

   std::string BuildOutputName(const std::string &OutputDir,
      const std::string &OutputPrefix,
      const std::string &ModeLabel,
      const std::string &IDLabel,
      const std::string &SideLabel,
      const std::string &OutputFormat)
   {
      return OutputDir + "/" + OutputPrefix + "_" + ModeLabel + "_" + IDLabel + "_" + SideLabel + "." + OutputFormat;
   }
}

int main(int argc, char *argv[])
{
   if(PrintHelpMessage(argc, argv))
      return 0;

   CommandLine CL(argc, argv);
   const std::string InputFileName = CL.Get("Input", "output/HFDistributions_First_100.root");
   const std::string OutputDir = CL.Get("OutputDir", "output");
   const std::string OutputPrefix = CL.Get("OutputPrefix", "HFDistributions");
   const std::string OutputFormat = CL.Get("OutputFormat", "pdf");
   const std::string HistogramSuffix = CL.Get("HistogramSuffix", "");
   const bool IgnoreEmptyWrapHigh = CL.GetBool("IgnoreEmptyWrapHigh", false);

   if(OutputFormat != "pdf" && OutputFormat != "png")
   {
      std::cerr << "Unsupported OutputFormat = " << OutputFormat << ". Allowed values are pdf or png." << std::endl;
      return 1;
   }

   TFile InputFile(InputFileName.c_str(), "READ");
   if(InputFile.IsZombie())
   {
      std::cerr << "Failed to open input file " << InputFileName << std::endl;
      return 1;
   }

   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gSystem->mkdir(OutputDir.c_str(), true);

   const std::vector<std::string> ModeLabels = {"CountCrossCheck", "AvgEnergy", "StdDevEnergy"};
   const std::vector<std::string> IDLabels = {"All", "id6", "id7"};
   bool FoundAtLeastOneSet = false;

   for(const std::string &ModeLabel : ModeLabels)
   {
      for(const std::string &IDLabel : IDLabels)
      {
         const std::string Suffix = "_" + ModeLabel + "_" + IDLabel;
         DistributionSet PositiveSet;
         DistributionSet NegativeSet;

         const bool HasPositive = LoadSet(InputFile, "Positive", Suffix, HistogramSuffix, PositiveSet, IgnoreEmptyWrapHigh);
         const bool HasNegative = LoadSet(InputFile, "Negative", Suffix, HistogramSuffix, NegativeSet, IgnoreEmptyWrapHigh);

         if(HasPositive == false && HasNegative == false)
         {
            std::cerr << "Missing histogram set for suffix " << Suffix << std::endl;
            return 1;
         }
         if(HasPositive == false || HasNegative == false)
         {
            std::cerr << "Incomplete histogram set for suffix " << Suffix << std::endl;
            DeleteSet(PositiveSet);
            DeleteSet(NegativeSet);
            return 1;
         }

         FoundAtLeastOneSet = true;
         const std::string PlusPDF = BuildOutputName(OutputDir, OutputPrefix, ModeLabel, IDLabel, "Plus", OutputFormat);
         const std::string MinusPDF = BuildOutputName(OutputDir, OutputPrefix, ModeLabel, IDLabel, "Minus", OutputFormat);

         DrawSidePDF(PlusPDF, "HF+ " + ModeLabel + " " + IDLabel, PositiveSet, ModeLabel);
         DrawSidePDF(MinusPDF, "HF- " + ModeLabel + " " + IDLabel, NegativeSet, ModeLabel);

         DeleteSet(PositiveSet);
         DeleteSet(NegativeSet);
      }
   }

   if(FoundAtLeastOneSet == false)
   {
      std::cerr << "Did not find any histograms in the expected new layout." << std::endl;
      return 1;
   }

   return 0;
}
