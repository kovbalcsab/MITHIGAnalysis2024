default: ExtractEmptyBXDistributons ExtractEmptyZBDistributons PrintInfoFromFile

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