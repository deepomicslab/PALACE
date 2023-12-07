# Compiler
CC = g++

# Compiler flags
CFLAGS = -Wall

# Libraries to link with
LDFLAGS = -L $(shell pwd)/libhts -lhts

all: $(shell pwd)/libhts/lib/libhts.a eref rDistance

$(shell pwd)/libhts/lib/libhts.a:
	cd htslib-1.18 && chmod u+x ./configure && ./configure prefix=$(shell pwd)/libhts && make && make install

eref: extract_ref.cpp $(shell pwd)/libhts/lib/libhts.a
	$(CC) $(CFLAGS) -o eref extract_ref.cpp $(LDFLAGS) -lpthread

rDistance: reads_overlap.cpp $(shell pwd)/libhts/lib/libhts.a
	$(CC) $(CFLAGS) -o rDistance reads_overlap.cpp $(LDFLAGS)

clean:
	cd htslib-1.18 && make clean
	rm -f eref rDistance