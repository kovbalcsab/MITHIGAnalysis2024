#ifndef INFOMANAGER_H
#define INFOMANAGER_H

#include <iostream>
#include <string>
#include <set>
#include <vector>

#include "TTimeStamp.h"
#include "TObject.h"
#include "TFile.h"
#include "TList.h"
#include "TMap.h"
#include "TParameter.h"
#include "TObjString.h"
#include "TString.h"
#include "TDirectory.h"

class GeneralInfoManager {
    public:
        GeneralInfoManager(TFile* File, TString Directory, bool ReadMode = false);
        ~GeneralInfoManager();

        void AddSourceFile(const std::string& filePath, const TTimeStamp* timeStamp);
        void AddCutParameter(const std::string& parameterName, double value, const TTimeStamp* timeStamp);

        TList* GetSourceFiles() const;
        std::vector<std::string> GetCreationTimes();
        TList* GetCutParameters() const;
        std::string GetSourceFileByTime(const std::string& timeStampStr) const;
        std::vector<std::string> GetParameterInfoFromString(const std::string& keyStr) const;
        double GetCutParameterValue(const std::string& parameterName, std::string timeStr = "");
        bool IsReadMode() const;
        bool IsInitialized() const;
        bool NoDuplicates() const;
        bool CutConsistencyCheck() const;
        bool NeedUpdatedLists() const;
        void UpdateLists();

        void ChangeReadMode(bool ReadMode);

        void SaveToFile();
        void PrintInfo() const;
    private:
        TFile* file;
        TString directory;
        TList* SourceFiles;
        std::vector<std::string> CreationTimes;
        TList* CutParameters;
        bool isReadMode;
        bool Initialized;
        bool needUpdatedLists;
};

#endif