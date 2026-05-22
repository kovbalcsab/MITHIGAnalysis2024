ProjectBase ?= ../..

default: ExecuteScanPFTowerLocations ExecuteMakeHFPFBinEdges ExecuteCheckHFPFBinning ExecutePlotHFDistributionsSides ExecuteHFDistributionExtraction ExecuteRecalculateMeanAndStd ExecuteDiffHFDistributions ExtractEmptyBXFullMaps PlotForestEtaPhiMaps

ExecuteScanPFTowerLocations: scan_pf_tower_locations.cpp
	g++ scan_pf_tower_locations.cpp -O3 -std=c++17 -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteScanPFTowerLocations `root-config --cflags --glibs`

ExecuteMakeHFPFBinEdges: MakeHFPFBinEdges.cpp
	g++ MakeHFPFBinEdges.cpp -O3 -std=c++17 -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteMakeHFPFBinEdges `root-config --cflags --glibs`

ExecuteCheckHFPFBinning: CheckHFPFBinning.cpp hf_pf_bin_edges.h
	g++ CheckHFPFBinning.cpp -O3 -std=c++17 -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteCheckHFPFBinning `root-config --cflags --glibs`

ExecutePlotHFDistributionsSides: PlotHFDistributionsSides.cpp
	g++ PlotHFDistributionsSides.cpp -O3 -std=c++17 -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecutePlotHFDistributionsSides `root-config --cflags --glibs`

ExecuteHFDistributionExtraction: manual/makeDistributions.cpp include/hfEdges_20260520.h include/utility.h
	g++ manual/makeDistributions.cpp -O3 -std=c++17 -I. -I$(ProjectBase)/CommonCode/include -I./include \
		-o ExecuteHFDistributionExtraction `root-config --cflags --glibs`

ExecuteRecalculateMeanAndStd: manual/recalculateMeanAndStd.cpp include/hfEdges_20260520.h include/utility.h
	g++ manual/recalculateMeanAndStd.cpp -O3 -std=c++17 -I. -I$(ProjectBase)/CommonCode/include -I./include \
		-o ExecuteRecalculateMeanAndStd `root-config --cflags --glibs`

ExecuteDiffHFDistributions: manual/diffHFDistributions.cpp include/hfEdges_20260520.h include/utility.h
	g++ manual/diffHFDistributions.cpp -O3 -std=c++17 -I. -I$(ProjectBase)/CommonCode/include -I./include \
		-o ExecuteDiffHFDistributions `root-config --cflags --glibs`

ExtractEmptyBXFullMaps: ForestEtaPhiMaps/ExtractEmptyBXFullMaps.C src/InfoManager.cpp
	g++ ForestEtaPhiMaps/ExtractEmptyBXFullMaps.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I./include \
		-o ExtractEmptyBXFullMaps `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

PlotForestEtaPhiMaps: ForestEtaPhiMaps/PlotForestEtaPhiMaps.C
	g++ ForestEtaPhiMaps/PlotForestEtaPhiMaps.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I./include \
		-o PlotForestEtaPhiMaps `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o
