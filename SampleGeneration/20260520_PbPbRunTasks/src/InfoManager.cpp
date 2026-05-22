#include "InfoManager.h"

GeneralInfoManager::GeneralInfoManager(TFile* File, TString Directory, bool ReadMode) {
    file = File;
    directory = Directory;
    isReadMode = ReadMode;
    Initialized = false;
    needUpdatedLists = true;
    CreationTimes = {};

    if (isReadMode) {
        SourceFiles = dynamic_cast<TList*>(file->Get(directory + "/SourceFiles"));
        CutParameters = dynamic_cast<TList*>(file->Get(directory + "/CutParameters"));

        if (SourceFiles && CutParameters) {
            Initialized = true;
        } else {
            std::cerr << "\033[1;31mError: Could not read all necessary information from the file.\033[0m" << std::endl;
        }

        UpdateLists();

        if (!NoDuplicates()) {
            std::cerr << "\033[1;33mWarning: Duplicates found in the information read from the file.\033[0m" << std::endl;
        }
        if (!CutConsistencyCheck()) {
            std::cerr << "\033[1;31mError: Inconsistent cut parameters found in the information read from the file.\033[0m" << std::endl;
        }
    } else {
        SourceFiles = new TList();
        SourceFiles->SetName("SourceFiles");
        SourceFiles->SetOwner(true);
        CutParameters = new TList();
        CutParameters->SetName("CutParameters");
        CutParameters->SetOwner(true);
        Initialized = true;
    }
}

GeneralInfoManager::~GeneralInfoManager() {
    if (!isReadMode) {
        delete SourceFiles;
        delete CutParameters;
    }
}

void GeneralInfoManager::AddSourceFile(const std::string& filePath, const TTimeStamp* timeStamp) {
    if (Initialized && !isReadMode) {
        SourceFiles->Add(new TObjString((std::string(timeStamp->AsString("s")) + "=" + filePath).c_str()));
    } else {
        std::cerr << "\033[1;31mError: Cannot add source file in read mode or if not initialized.\033[0m" << std::endl;
    }

    needUpdatedLists = true;
}

void GeneralInfoManager::AddCutParameter(const std::string& parameterName, double value, const TTimeStamp* timeStamp) {
    if (Initialized && !isReadMode) {
        CutParameters->Add(new TObjString((parameterName+"("+std::string(timeStamp->AsString("s"))+")="+std::to_string(value)).c_str()));
    } else {
        std::cerr << "\033[1;31mError: Cannot add cut parameter in read mode or if not initialized.\033[0m" << std::endl;
    }
}

TList* GeneralInfoManager::GetSourceFiles() const {
    return SourceFiles;
}

std::vector<std::string> GeneralInfoManager::GetCreationTimes() {
    if (needUpdatedLists) {
        UpdateLists();
    }
    return CreationTimes;
}

TList* GeneralInfoManager::GetCutParameters() const {
    return CutParameters;
}

bool GeneralInfoManager::NeedUpdatedLists() const {
    return needUpdatedLists;
}

void GeneralInfoManager::UpdateLists() {
    if (!Initialized) {
        std::cerr << "\033[1;31mError: Cannot update lists if not initialized.\033[0m" << std::endl;
        return;
    }

    // Update CreationTimes list based on the keys in SourceFiles
    CreationTimes.clear();

    TIter next(SourceFiles);
    TObjString* current;
    while ((current = dynamic_cast<TObjString*>(next()))) {
        CreationTimes.push_back(std::string(current->GetString().Data()).substr(0, std::string(current->GetString().Data()).find('=')));
    }

    needUpdatedLists = false;
}

std::string GeneralInfoManager::GetSourceFileByTime(const std::string& timeStampStr) const {
    if (!Initialized) {
        std::cerr << "\033[1;31mError: Cannot get source file by time if not initialized.\033[0m" << std::endl;
        return ""; // Return a default value or consider throwing an exception
    }

    TIter next(SourceFiles);
    TObjString* current;
    while ((current = dynamic_cast<TObjString*>(next()))) {
        if (std::string(current->GetString().Data()).substr(0, std::string(current->GetString().Data()).find('=')) == timeStampStr) {
            return std::string(current->GetString().Data()).substr(std::string(current->GetString().Data()).find('=') + 1);
        }
    }

    
    std::cerr << "\033[1;31mError: Source file not found for the given time stamp.\033[0m" << std::endl;
    return ""; // Return a default value or consider throwing an exception
}

std::vector<std::string> GeneralInfoManager::GetParameterInfoFromString(const std::string& keyStr) const {
    std::vector<std::string> result;
    if (!Initialized) {
        std::cerr << "\033[1;31mError: Cannot get parameter info from string if not initialized.\033[0m" << std::endl;
        return result; // Return an empty vector or consider throwing an exception
    }

    result.push_back(keyStr.substr(0, keyStr.find('('))); // Add the parameter name without the timestamp
    result.push_back(keyStr.substr(keyStr.find('(') + 1, keyStr.find(')') - keyStr.find('(') - 1)); // Add the timestamp
    result.push_back(keyStr.substr(keyStr.find('=') + 1)); // Add the parameter value
    return result;
}

double GeneralInfoManager::GetCutParameterValue(const std::string& parameterName, std::string timeStampStr) {
    if (!Initialized) {
        std::cerr << "\033[1;31mError: Cannot get cut parameter value if not initialized.\033[0m" << std::endl;
        return -977; // Return a default value or consider throwing an exception
    }

    if (needUpdatedLists) {
        UpdateLists();
    }

    if (timeStampStr == "") {
        timeStampStr = CreationTimes.size() > 0 ? CreationTimes[0] : "";
    }
    TIter next(CutParameters);
    TObjString* current;
    while ((current = dynamic_cast<TObjString*>(next()))) {
        std::string currentKey = current->GetString().Data();
        std::vector<std::string> paramInfo = GetParameterInfoFromString(currentKey);
        if (paramInfo.size() == 3 && paramInfo[0] == parameterName && paramInfo[1] == timeStampStr) {
            try {
                return std::stod(paramInfo[2]);
            } catch (const std::invalid_argument& e) {
                std::cerr << "\033[1;31mError: Invalid parameter value format for '" << currentKey << "'.\033[0m" << std::endl;
                return -977; // Return a default value or consider throwing an exception
            }
        }
    }

    std::cerr << "\033[1;31mError: Cut parameter '" << parameterName << "(" << timeStampStr << ")' not found.\033[0m" << std::endl;
    return -977; // Return a default value or consider throwing an exception
}

bool GeneralInfoManager::IsReadMode() const {
    return isReadMode;
}

bool GeneralInfoManager::IsInitialized() const {
    return Initialized;
}

bool GeneralInfoManager::NoDuplicates() const {
    if (!Initialized) {
        std::cerr << "\033[1;31mError: Cannot check for duplicates if not initialized.\033[0m" << std::endl;
        return false;
    }

    // Check for duplicates in SourceFiles
    TIter next(SourceFiles);
    TObjString* current;
    std::set<std::string> uniqueSourceFiles;
    while ((current = dynamic_cast<TObjString*>(next()))) {
        if (uniqueSourceFiles.find(current->GetString().Data()) != uniqueSourceFiles.end()) {
            return false; // Duplicate found
        }
        uniqueSourceFiles.insert(current->GetString().Data());
    }
    return true; // No duplicates found  
}

bool GeneralInfoManager::CutConsistencyCheck() const {
    if (!Initialized) {
        std::cerr << "\033[1;31mError: Cannot check cut consistency if not initialized.\033[0m" << std::endl;
        return false;
    }

    //Check if a cut parameter has only one value, possibly multiple occurances, but no conflicting values.
    TIter next(CutParameters);
    TObjString* current;
    std::map<std::string, std::set<std::string>> parameterValues; // Map to store parameter names and their associated values
    while ((current = dynamic_cast<TObjString*>(next()))) {
        std::vector<std::string> paramInfo = GetParameterInfoFromString(current->GetString().Data());
        if (paramInfo.size() != 3) {
            std::cerr << "\033[1;31mError: Invalid cut parameter format for '" << current->GetString().Data() << "'.\033[0m" << std::endl;
            continue; // Skip this entry and continue checking others
        }
        if (paramInfo.size() == 3) {
            // Check if the parameter name already exists in the map
            if (parameterValues.find(paramInfo[0]) != parameterValues.end()) {
                // If it exists, check if the value is different from the existing values
                if (parameterValues[paramInfo[0]].find(paramInfo[2]) == parameterValues[paramInfo[0]].end()) {
                    std::cerr << "\033[1;31mError: Inconsistent values found for cut parameter '" << paramInfo[0] << "'.\033[0m" << std::endl;
                    return false; // Inconsistent values found
                }
            }
            std::string paramName = paramInfo[0];
            std::string paramValue = paramInfo[2];
            parameterValues[paramName].insert(paramValue);
        }
    }

    return true; // Assuming all checks passed for now
}

void GeneralInfoManager::ChangeReadMode(bool ReadMode) {
    isReadMode = ReadMode;
    std::cerr << "\033[1;33mWarning: Changing read mode after initialization may lead to inconsistent state.\033[0m" << std::endl;
}

void GeneralInfoManager::SaveToFile() {
    if (Initialized && !isReadMode) {
        file->cd();
        TDirectory* dir = file->GetDirectory(directory);
        if (!dir) {
            dir = file->mkdir(directory);
        }
        dir->cd();
        dir->WriteTObject(SourceFiles);
        dir->WriteTObject(CutParameters);

        file->Write();
    } else {
        std::cerr << "\033[1;31mError: Cannot save to file in read mode or if not initialized.\033[0m" << std::endl;
    }
}

void GeneralInfoManager::PrintInfo() const {
    if (!Initialized) {
        std::cerr << "\033[1;31mError: Cannot print info if not initialized.\033[0m" << std::endl;
        return;
    }

    std::cout << "Source Files:" << std::endl;
    TIter next(SourceFiles);
    TObjString* current;
    while ((current = dynamic_cast<TObjString*>(next()))) {
        std::cout << " - " << (std::string(current->GetString().Data()).substr(0, std::string(current->GetString().Data()).find('='))) << ": " << std::string(current->GetString().Data()).substr(std::string(current->GetString().Data()).find('=') + 1) << std::endl;
    }
    int nSourceFiles = SourceFiles->GetEntries();
    std::cout << "Total Source Files: " << nSourceFiles << std::endl;

    std::cout << "Cut Parameters:" << std::endl;
    TIter paramNext(CutParameters);
    TObjString* currentParam;
    while ((currentParam = dynamic_cast<TObjString*>(paramNext()))) {
        std::vector<std::string> paramInfo = GetParameterInfoFromString(currentParam->GetString().Data());
        if (paramInfo.size() == 3) {
            if (nSourceFiles > 1) {
                std::cout << " - " << paramInfo[0] << " (" << paramInfo[1] << "): " << paramInfo[2] << std::endl;
            } else {
                std::cout << " - " << paramInfo[0] << ": " << paramInfo[2] << std::endl;
            }
        } else {
            std::cout << "\033[1;31m - " << currentParam->GetString().Data() << "\033[0m" << std::endl;
        }
    }
}
