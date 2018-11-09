prefix=$(shell pwd)

BAMTOOLS_INCLUDE_DIR=$(prefix)/bamtools/install/include/bamtools
BAMTOOLS_LIB_DIR=$(prefix)/bamtools/install/lib/bamtools

all:
	g++ -std=c++14 -o Blacklist blacklist.cpp -I$(BAMTOOLS_INCLUDE_DIR) -L$(BAMTOOLS_LIB_DIR) -lbamtools -lz -Wl,-rpath,$(BAMTOOLS_LIB_DIR)

debug:
	g++ -std=c++14 -g -o Blacklist blacklist.cpp -I$(BAMTOOLS_INCLUDE_DIR) -L$(BAMTOOLS_LIB_DIR) -lbamtools -lz -Wl,-rpath,$(BAMTOOLS_LIB_DIR)

blacklist:
	./Blacklist chr1 > final.chr1.out
	./Blacklist chr2 > final.chr2.out
	./Blacklist chr3 > final.chr3.out
	./Blacklist chr4 > final.chr4.out
	./Blacklist chr5 > final.chr5.out
	./Blacklist chr6 > final.chr6.out
	./Blacklist chr7 > final.chr7.out
	./Blacklist chr8 > final.chr8.out
	./Blacklist chr9 > final.chr9.out
	./Blacklist chr10 > final.chr10.out
	./Blacklist chr11 > final.chr11.out
	./Blacklist chr12 > final.chr12.out
	./Blacklist chr13 > final.chr13.out
	./Blacklist chr14 > final.chr14.out
	./Blacklist chr15 > final.chr15.out
	./Blacklist chr16 > final.chr16.out
	./Blacklist chr17 > final.chr17.out
	./Blacklist chr18 > final.chr18.out
	./Blacklist chr19 > final.chr19.out
	./Blacklist chr20 > final.chr20.out
	./Blacklist chr21 > final.chr21.out
	./Blacklist chr22 > final.chr22.out
	./Blacklist chrX > final.chrX.out
	./Blacklist chrY > final.chrY.out

