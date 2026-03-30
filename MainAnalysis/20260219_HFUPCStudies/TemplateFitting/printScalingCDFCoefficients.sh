#!/usr/bin/env bash
set -euo pipefail

if [ "${1:-}" = "" ]; then
  echo "Usage: $0 <ScalingCDFResult.root> [FunctionName]" >&2
  exit 1
fi

RootFile="$1"
FunctionName="${2:-fitFunc}"

if [ ! -f "$RootFile" ]; then
  echo "Error: ROOT file not found: $RootFile" >&2
  exit 1
fi

ROOT_FILE="$RootFile" FUNCTION_NAME="$FunctionName" root -l -b <<'ROOT_MACRO'
{
  const char *fileName = gSystem->Getenv("ROOT_FILE");
  const char *functionName = gSystem->Getenv("FUNCTION_NAME");

  TFile *file = TFile::Open(fileName, "READ");
  if(file == nullptr || file->IsZombie())
  {
    std::cerr << "Error: could not open file " << fileName << std::endl;
    gSystem->Exit(2);
  }

  TF1 *func = dynamic_cast<TF1 *>(file->Get(functionName));
  if(func == nullptr)
  {
    std::cerr << "Error: function '" << functionName << "' was not found in " << fileName << std::endl;
    file->Close();
    gSystem->Exit(3);
  }

  std::cout << "File: " << fileName << std::endl;
  std::cout << "Function: " << func->GetName() << std::endl;
  std::cout << "Formula: " << func->GetExpFormula("p").Data() << std::endl;
  if(func->GetNDF() > 0)
    std::cout << "Chi2/NDF: " << func->GetChisquare() << " / " << func->GetNDF()
              << " = " << (func->GetChisquare() / func->GetNDF()) << std::endl;
  std::cout << "Coefficients:" << std::endl;
  for(int i = 0; i < func->GetNpar(); ++i)
  {
    const char *parName = func->GetParName(i);
    if(parName == nullptr || parName[0] == '\0')
      parName = Form("p%d", i);
    std::cout << "  " << parName << " = " << func->GetParameter(i)
              << " +/- " << func->GetParError(i) << std::endl;
  }

  file->Close();
  gSystem->Exit(0);
}
ROOT_MACRO
