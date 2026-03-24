#ifndef ROOT_IO_UTILS_H
#define ROOT_IO_UTILS_H

#include <iostream>
#include <string>

#include "TBranch.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TTree.h"

namespace RootIOUtils
{
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
    if(file != nullptr)
    {
      file->Close();
      delete file;
    }
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
} // namespace RootIOUtils

#endif

