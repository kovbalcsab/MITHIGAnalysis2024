default: ExtractEmptyBXDistributons ExtractEmptyZBDistributons PrintInfoFromFile ExtractSignalMCDistributons PlotComparisonMultiVar PlotForestEtaPhiMaps ScalingFit ScalingCDF plotScalingResult PlotCDFResult ExtractEmptyBXFullMaps

ExtractEmptyBXDistributons: EmptyBX/ExtractEmptyBXDistributons.C src/InfoManager.cpp
	g++ EmptyBX/ExtractEmptyBXDistributons.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ExtractEmptyBXDistributons `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

ExtractEmptyZBDistributons: EmptyZB/ExtractEmptyZBDistributons.C src/InfoManager.cpp
	g++ EmptyZB/ExtractEmptyZBDistributons.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ExtractEmptyZBDistributons `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

PrintInfoFromFile: printInfoFromFile.C src/InfoManager.cpp
	g++ printInfoFromFile.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o PrintInfoFromFile `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

ExtractSignalMCDistributons: MCSignal/ExtractSignalMCDistributons.C src/InfoManager.cpp
	g++ MCSignal/ExtractSignalMCDistributons.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ExtractSignalMCDistributons `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

PlotComparisonMultiVar: Plotting/PlotComparisonMultiVar.C
	g++ Plotting/PlotComparisonMultiVar.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o PlotComparisonMultiVar `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

PlotForestEtaPhiMaps: Plotting/PlotForestEtaPhiMaps.C
	g++ Plotting/PlotForestEtaPhiMaps.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o PlotForestEtaPhiMaps `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

ScalingFit: TemplateFitting/ScalingFit.C
	g++ TemplateFitting/ScalingFit.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ScalingFit `root-config --cflags --glibs` -lRooFit -lRooFitCore \
		$(ProjectBase)/CommonCode/library/Messenger.o

ScalingCDF: TemplateFitting/ScalingCDF.C
	g++ TemplateFitting/ScalingCDF.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ScalingCDF `root-config --cflags --glibs` -lRooFit -lRooFitCore \
		$(ProjectBase)/CommonCode/library/Messenger.o

plotScalingResult: TemplateFitting/plotScalingResult.C
	g++ TemplateFitting/plotScalingResult.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o plotScalingResult `root-config --cflags --glibs` -lRooFit -lRooFitCore \
		$(ProjectBase)/CommonCode/library/Messenger.o

PlotCDFResult: TemplateFitting/plotCDFResult.C
	g++ TemplateFitting/plotCDFResult.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o PlotCDFResult `root-config --cflags --glibs` -lRooFit -lRooFitCore \
		$(ProjectBase)/CommonCode/library/Messenger.o

ExtractEmptyBXFullMaps: ForestEtaPhiMaps/ExtractEmptyBXFullMaps.C src/InfoManager.cpp
	g++ ForestEtaPhiMaps/ExtractEmptyBXFullMaps.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ExtractEmptyBXFullMaps `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o