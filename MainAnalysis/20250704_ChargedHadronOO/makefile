default: ExecuteChargedHadronRAA ExecuteSystematicPlotter

ExecuteChargedHadronRAA: ChargedHadronRAA.cpp
	g++ ChargedHadronRAA.cpp -O3 -I. -I$(ProjectBase)/CommonCode/include \
	  -I./include \
		-o ExecuteChargedHadronRAA `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o

ExecuteSystematicPlotter: SystematicPlotter.cpp
	g++ SystematicPlotter.cpp -O3 -I. -I$(ProjectBase)/CommonCode/include \
	  -I./include \
		-o ExecuteSystematicPlotter `root-config --cflags --glibs` \
		$(ProjectBase)/CommonCode/library/Messenger.o
		