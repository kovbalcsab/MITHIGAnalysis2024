#ifndef BX_UTILS
#define BX_UTILS

#include <fstream>

enum bxSchemeBits { kEmpty = 0, kBx1 = 1, kBx2 = 2, kBxBoth = 3 };

map<int, bxSchemeBits> getBxScheme(string schemeCSV) {
  map<int, bxSchemeBits> bxScheme;

  ifstream inFile(schemeCSV.c_str());
  string line;
  while (getline(inFile, line)) {
    // std::cout << "LINE: " << line << std::endl;

    if (line.size() == 0)
      continue;
    if (line.find("bunch_number") != std::string::npos)
      continue;
    // Processing csv so comma separate
    vector<string> lineCommaSep;
    while (line.find(",") != string::npos && line.find(",") != 0) {
      lineCommaSep.push_back(line.substr(0, line.find(",")));
      line.replace(0, line.find(",") + 1, "");
    }

    // Grab bxnum, and define the bit as empty, bx1, bx2, or both
    int bxNum = std::stoi(lineCommaSep[0]);
    double bx1 = std::stod(lineCommaSep[2]);
    double bx2 = std::stod(lineCommaSep[3]);
    bxSchemeBits bxBit = kEmpty;
    if (bx1 > 0.0 && bx2 > 0.0)
      bxBit = kBxBoth;
    else if (bx1 > 0.0)
      bxBit = kBx1;
    else if (bx2 > 0.0)
      bxBit = kBx2;

    bxScheme[bxNum] = bxBit;
  }
  inFile.close();

  return bxScheme;
}

bool lookForBX(int bx_num, map<int, bxSchemeBits> &bxScheme, int BX_Sel) {
  if (BX_Sel == 0)
    return true; // No BX selection
  else if (BX_Sel > 0) {
    for (int i = 1; i <= BX_Sel; i++) {
      if (bxScheme[bx_num - i] != kEmpty)
        return false; // Not an empty BX before
    }
    return true; // Empty BX before
  } else {
    for (int i = 1; i <= -BX_Sel; i++) {
      if (bxScheme[bx_num + i] != kEmpty)
        return false; // Not an empty BX after
    }
    return true; // Empty BX after
  }
}

#endif // BX_UTILS
