ProjectBase ?= ../..

default: ExecuteScanPFTowerLocations ExecuteMakeHFPFBinEdges ExecuteCheckHFPFBinning ExecutePlotHFDistributionsSides ExecuteHFDistributionExtraction

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

ExecuteHFDistributionExtraction: manual/makeDistributions.cpp manual/hfEdges_20260520.h
	g++ manual/makeDistributions.cpp -O3 -std=c++17 -I. -I$(ProjectBase)/CommonCode/include \
		-o ExecuteHFDistributionExtraction `root-config --cflags --glibs`
