#ifndef PODIFFRACTIVESELECTIONFUNCS
#define PODIFFRACTIVESELECTIONFUNCS

bool getBaselinepOEventSel(pODiffractiveTreeMessenger &ch) {
    // zVtx cut
    if (fabs(ch.VZ) >= 15.){ return false; }

    // Cluster Compatibility Filter
    if (ch.ClusterCompatibilityFilter == false) { return false; }

    // PV selection
    if (ch.PVFilter == false) { return false; }

    return true;
}

#endif
