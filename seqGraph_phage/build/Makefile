# Compiler
CC = g++

# Compiler flags
CFLAGS = -Wall -Ihtslib-1.18

# Libraries to link with
LDFLAGS = -Lhtslib-1.18 -lhts

all: htslib-1.18/libhts.a eref rDistance

htslib-1.18/libhts.a:
	cd htslib-1.18 && chmod u+x ./configure && ./configure && make

eref: extract_ref.cpp htslib-1.18/libhts.a
	$(CC) $(CFLAGS) -o eref extract_ref.cpp $(LDFLAGS)

rDistance: reads_overlap.cpp reads_overlap.h htslib-1.18/libhts.a
	$(CC) $(CFLAGS) -o rDistance reads_overlap.cpp $(LDFLAGS)

clean:
	cd htslib-1.18 && make clean
	rm -f eref rDistance
