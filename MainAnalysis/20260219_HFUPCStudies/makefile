default: ExtractEmptyBXDistributons ExtractEmptyMCDistributons ExtractEmptyZBDistributons PrintInfoFromFile ExtractSignalMCDistributons PlotComparisonMultiVar PlotComparisonEnergyDistributions PlotForestEtaPhiMaps ScalingFit ScalingCDF plotScalingResult plotScalingResultHist PlotCDFResult ExtractEmptyBXFullMaps ExtractUPCDataSignalDistributons TemplateFitting/Unfold TemplateFitting/Unfold_SVD TemplateFitting/unfolding_iter_optimization TemplateFitting/Deconvolve_FFT TemplateFitting/deconvolve_fft_validation TemplateFitting/unfold_svd_validation TemplateFitting/plot_validation_method_comparison TemplateFitting/plot_validation_quarter_chi2

ExtractEmptyBXDistributons: EmptyBX/ExtractEmptyBXDistributons.C src/InfoManager.cpp
	g++ EmptyBX/ExtractEmptyBXDistributons.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ExtractEmptyBXDistributons `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

ExtractEmptyMCDistributons: EmptyMC/ExtractEmptyMCDistributons.C src/InfoManager.cpp
	g++ EmptyMC/ExtractEmptyMCDistributons.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ExtractEmptyMCDistributons `root-config --cflags --glibs` \
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

PlotComparisonEnergyDistributions: Plotting/PlotComparisonEnergyDistributions.C
	g++ Plotting/PlotComparisonEnergyDistributions.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o PlotComparisonEnergyDistributions `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

PlotForestEtaPhiMaps: Plotting/PlotForestEtaPhiMaps.C
	g++ Plotting/PlotForestEtaPhiMaps.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o PlotForestEtaPhiMaps `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

ScalingFit: TemplateFitting/ScalingFit.C src/InfoManager.cpp
	g++ TemplateFitting/ScalingFit.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ScalingFit `root-config --cflags --glibs` -lRooFit -lRooFitCore \
		$(ProjectBase)/CommonCode/library/Messenger.o

ScalingCDF: TemplateFitting/ScalingCDF.C src/InfoManager.cpp
	g++ TemplateFitting/ScalingCDF.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
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

plotScalingResultHist: TemplateFitting/plotScalingResultHist.C
	g++ TemplateFitting/plotScalingResultHist.C -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o plotScalingResultHist `root-config --cflags --glibs` -lRooFit -lRooFitCore \
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

ExtractUPCDataSignalDistributons: UPCDataSignal/ExtractUPCDataSignalDistributons.C src/InfoManager.cpp
	g++ UPCDataSignal/ExtractUPCDataSignalDistributons.C src/InfoManager.cpp -I. -O3 -I$(ProjectBase)/CommonCode/include \
	  -I$(ProjectBase)/WeightHandler/20250305_DzeroUPC_GptGyWeight/include \
	  -I$(ProjectBase)/WeightHandler/20250312_DzeroUPC_multiplicityWeight/include \
	  -I./include \
		-o ExtractUPCDataSignalDistributons `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

TemplateFitting/Unfold: TemplateFitting/Unfold.C src/InfoManager.cpp
	g++ TemplateFitting/Unfold.C src/InfoManager.cpp -o TemplateFitting/Unfold `root-config --cflags --glibs` -O3 -I. -I../../CommonCode/include \
		-I$(HOME)/RooUnfold/src -L$(HOME)/RooUnfold -lRooUnfold -I./include

TemplateFitting/Unfold_SVD: TemplateFitting/Unfold_SVD.C src/InfoManager.cpp
	g++ TemplateFitting/Unfold_SVD.C src/InfoManager.cpp -o TemplateFitting/Unfold_SVD `root-config --cflags --glibs` -O3 -I. -I../../CommonCode/include \
		-I$(HOME)/RooUnfold/src -L$(HOME)/RooUnfold -lRooUnfold -I./include

TemplateFitting/Deconvolve_FFT: TemplateFitting/Deconvolve_FFT.C src/InfoManager.cpp
	g++ TemplateFitting/Deconvolve_FFT.C src/InfoManager.cpp -o TemplateFitting/Deconvolve_FFT `root-config --cflags --glibs` -O3 -I. -I../../CommonCode/include \
		-I$(HOME)/RooUnfold/src -L$(HOME)/RooUnfold -lRooUnfold -I./include

TemplateFitting/unfolding_iter_optimization: TemplateFitting/unfolding_iter_optimization.C src/InfoManager.cpp
	g++ TemplateFitting/unfolding_iter_optimization.C src/InfoManager.cpp -o TemplateFitting/unfolding_iter_optimization \
		`root-config --cflags --glibs` -O3 -I. -I../../CommonCode/include -I$(HOME)/RooUnfold/src -L$(HOME)/RooUnfold -lRooUnfold -I./include

TemplateFitting/deconvolve_fft_validation: TemplateFitting/deconvolve_fft_validation.C src/InfoManager.cpp
	g++ TemplateFitting/deconvolve_fft_validation.C src/InfoManager.cpp -o TemplateFitting/deconvolve_fft_validation \
		`root-config --cflags --glibs` -O3 -I. -I../../CommonCode/include -I$(HOME)/RooUnfold/src -L$(HOME)/RooUnfold -lRooUnfold -I./include

TemplateFitting/unfold_svd_validation: TemplateFitting/unfold_svd_validation.C src/InfoManager.cpp
	g++ TemplateFitting/unfold_svd_validation.C src/InfoManager.cpp -o TemplateFitting/unfold_svd_validation \
		`root-config --cflags --glibs` -O3 -I. -I../../CommonCode/include -I$(HOME)/RooUnfold/src -L$(HOME)/RooUnfold -lRooUnfold -I./include

TemplateFitting/plot_validation_method_comparison: TemplateFitting/plot_validation_method_comparison.C
	g++ TemplateFitting/plot_validation_method_comparison.C -o TemplateFitting/plot_validation_method_comparison \
		`root-config --cflags --glibs` -O3 -I. -I../../CommonCode/include -I./include

TemplateFitting/plot_validation_quarter_chi2: TemplateFitting/plot_validation_quarter_chi2.C
	g++ TemplateFitting/plot_validation_quarter_chi2.C -o TemplateFitting/plot_validation_quarter_chi2 \
		`root-config --cflags --glibs` -O3 -I. -I../../CommonCode/include -I./include
