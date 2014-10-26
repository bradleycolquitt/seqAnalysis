#Makefile


CC=g++
BAMTOOLS_ROOT=/usr/local/bamtools
CXXFLAGS=-v -Wl,-rpath,$(BAMTOOLS_ROOT)/lib

test1: 
	$(CC) $(CXXFLAGS) -c test.cpp -o bamtools_test -I /usr/local/bamtools/include -L /usr/local/bamtools/lib -lbamtools -lz 

test2: 
	$(CC) $(CXXFLAGS) -c test2.cpp  -I /usr/local/bamtools/include -L /usr/local/bamtools/lib -lbamtools -lz -o bamtools_test2
