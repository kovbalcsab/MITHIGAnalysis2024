#ifndef ROOT_IO_UTILS_H
#define ROOT_IO_UTILS_H

#include <cerrno>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <string>

#include "TBranch.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLeaf.h"
#include "TList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TObject.h"
#include "TTree.h"

namespace RootIOUtils
{
inline void CloseAndDeleteFile(TFile *file)
{
  if(file == nullptr)
    return;
  if(file->IsOpen())
    file->Close();
  delete file;
}

inline void PrintAvailableBranches(TTree *tree)
{
  if(tree == nullptr)
  {
    std::cerr << "(tree is null)" << std::endl;
    return;
  }

  TObjArray *branches = tree->GetListOfBranches();
  if(branches == nullptr || branches->GetEntries() == 0)
  {
    std::cerr << "(no branches)" << std::endl;
    return;
  }

  for(int i = 0; i < branches->GetEntries(); ++i)
  {
    TObject *obj = branches->At(i);
    if(obj != nullptr)
      std::cerr << obj->GetName() << (i + 1 < branches->GetEntries() ? ", " : "\n");
  }
}

inline TFile *OpenFileOrNull(const std::string &fileName, const char *mode, const std::string &label)
{
  TFile *file = TFile::Open(fileName.c_str(), mode);
  if(file == nullptr || file->IsZombie())
  {
    std::cerr << "Error opening " << label << ": " << fileName << std::endl;
    CloseAndDeleteFile(file);
    return nullptr;
  }
  return file;
}

inline TTree *GetTreeOrNull(TFile *file, const std::string &treeName, const std::string &label)
{
  if(file == nullptr)
  {
    std::cerr << "Internal error: null file while loading tree '" << treeName << "' from " << label << std::endl;
    return nullptr;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get(treeName.c_str()));
  if(tree == nullptr)
    std::cerr << "Could not find tree '" << treeName << "' in " << label << std::endl;
  return tree;
}

inline TBranch *RequireBranchOrNull(TTree *tree, const std::string &branchName, const std::string &label)
{
  if(tree == nullptr)
  {
    std::cerr << "Internal error: null tree while looking for branch '" << branchName << "' in " << label << std::endl;
    return nullptr;
  }

  TBranch *branch = tree->GetBranch(branchName.c_str());
  if(branch == nullptr)
  {
    std::cerr << "Could not find branch: " << branchName << " in " << label << std::endl;
    std::cerr << "Available branches: ";
    PrintAvailableBranches(tree);
  }
  return branch;
}

template <typename TObjectType>
inline TObjectType *GetObjectOrNull(TFile *file, const std::string &objectName, const std::string &label)
{
  if(file == nullptr)
  {
    std::cerr << "Internal error: null file while loading object '" << objectName << "' from " << label << std::endl;
    return nullptr;
  }

  TObjectType *object = dynamic_cast<TObjectType *>(file->Get(objectName.c_str()));
  if(object == nullptr)
    std::cerr << "Could not find object '" << objectName << "' in " << label << std::endl;
  return object;
}

inline TH1D *LoadTreeBranchHistogramOrNull(const std::string &fileName,
                                           const std::string &treeName,
                                           const std::string &branchName,
                                           const std::string &histName,
                                           const std::string &histTitle,
                                           int nBins, double xMin, double xMax,
                                           int dataQuarter = -1,
                                           bool printLoadSummary = true,
                                           bool useSumw2 = false,
                                           const std::string &label = "input file")
{
  TFile *file = OpenFileOrNull(fileName, "READ", label);
  if(file == nullptr)
    return nullptr;

  TTree *tree = GetTreeOrNull(file, treeName, label + ": " + fileName);
  if(tree == nullptr)
  {
    CloseAndDeleteFile(file);
    return nullptr;
  }

  TBranch *branch = RequireBranchOrNull(tree, branchName, label + ": " + treeName);
  if(branch == nullptr)
  {
    CloseAndDeleteFile(file);
    return nullptr;
  }

  TLeaf *leaf = branch->GetLeaf(branchName.c_str());
  if(leaf == nullptr)
    leaf = static_cast<TLeaf *>(branch->GetListOfLeaves()->First());
  const std::string leafType = (leaf != nullptr && leaf->GetTypeName() != nullptr) ? leaf->GetTypeName() : "";
  const bool isFloat = (leafType == "Float_t" || leafType == "float");

  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);
  TH1D *hist = new TH1D(histName.c_str(), histTitle.c_str(), nBins, xMin, xMax);
  TH1::AddDirectory(oldAddDirectory);
  if(useSumw2)
    hist->Sumw2();

  double valueDouble = 0.0;
  float valueFloat = 0.0f;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus(branchName.c_str(), 1);
  if(isFloat)
    tree->SetBranchAddress(branchName.c_str(), &valueFloat);
  else
    tree->SetBranchAddress(branchName.c_str(), &valueDouble);

  const Long64_t nEntries = tree->GetEntries();
  Long64_t nSelected = 0;
  for(Long64_t i = 0; i < nEntries; ++i)
  {
    if(dataQuarter >= 0 && (i % 4) != dataQuarter)
      continue;
    tree->GetEntry(i);
    hist->Fill(isFloat ? static_cast<double>(valueFloat) : valueDouble);
    nSelected++;
  }

  tree->SetBranchStatus("*", 1);
  tree->ResetBranchAddresses();
  CloseAndDeleteFile(file);

  if(printLoadSummary)
  {
    if(dataQuarter >= 0)
      std::cout << "  Loaded " << nSelected << "/" << nEntries << " entries from "
                << fileName << "/" << treeName << "/" << branchName
                << " (quarter " << dataQuarter << "/4)" << std::endl;
    else
      std::cout << "  Loaded " << nSelected << " entries from "
                << fileName << "/" << treeName << "/" << branchName << std::endl;
  }
  return hist;
}

inline int InferCutParameterIntOrDefault(TFile *file, const std::string &parameterName,
                                         int defaultValue = -1,
                                         const std::string &cutListPath = "InfoDir/CutParameters")
{
  if(file == nullptr)
    return defaultValue;

  TList *cuts = dynamic_cast<TList *>(file->Get(cutListPath.c_str()));
  if(cuts == nullptr)
    return defaultValue;

  const std::string prefix = parameterName + "(";
  TIter next(cuts);
  TObject *obj = nullptr;
  while((obj = next()) != nullptr)
  {
    TObjString *entry = dynamic_cast<TObjString *>(obj);
    if(entry == nullptr)
      continue;

    const std::string text = entry->GetString().Data();
    if(text.rfind(prefix, 0) != 0)
      continue;

    const std::size_t eqPos = text.find('=');
    if(eqPos == std::string::npos)
      return defaultValue;

    const char *valueText = text.c_str() + eqPos + 1;
    char *endPtr = nullptr;
    errno = 0;
    const double parsedValue = std::strtod(valueText, &endPtr);
    while(endPtr != nullptr && std::isspace(static_cast<unsigned char>(*endPtr)))
      ++endPtr;
    if(valueText == endPtr || errno != 0 || (endPtr != nullptr && *endPtr != '\0'))
      return defaultValue;

    return static_cast<int>(std::lround(parsedValue));
  }

  return defaultValue;
}
} // namespace RootIOUtils

#endif
