Default: Execute Shuffle

Execute: ReduceForest.cpp
	g++ ReduceForest.cpp -o Execute \
		`root-config --cflags --glibs` \
		-I$(ProjectBase)/CommonCode/include $(ProjectBase)/CommonCode/library/Messenger.o

Shuffle: ShuffleSkims.cpp
	g++ ShuffleSkims.cpp -o Shuffle \
		`root-config --cflags --glibs` \
		-I$(ProjectBase)/CommonCode/include $(ProjectBase)/CommonCode/library/Messenger.o
