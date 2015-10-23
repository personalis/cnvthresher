# Makefile for cnvthresher
# You MUST change BAMTOOLS to point to your bamtools installation
BAMTOOLS=/PATH/TO/YOUR/bamtools-2.3.0
CC=g++
INCLUDES=-I$(BAMTOOLS)/include
CLIBS=-L$(BAMTOOLS)/lib -l:libbamtools.so.2.3.0 
CFLAGS=-c -g -Wall -ansi -D_GNU_SOURCE
OBJECTS=cnvthresher.o GFF.o GRegion.o FileUtil.o BamEvidence.o EdgeEvidence.o InsSizeEvidence.o MorphEvidence.o Fasta.o
all: cnvthresher

cnvthresher: $(OBJECTS)
	$(CC) $(OBJECTS) -o cnvthresher $(CLIBS)

cnvthresher.o: cnvthresher.cpp
	$(CC) $(INCLUDES) $(CFLAGS) cnvthresher.cpp

GFF.o: GFF.cpp GFF.h
	$(CC) $(INCLUDES) $(CFLAGS) GFF.cpp

GRegion.o: GRegion.cpp GRegion.h
	$(CC) $(INCLUDES) $(CFLAGS) GRegion.cpp

FileUtil.o: FileUtil.cpp FileUtil.h
	$(CC) $(INCLUDES) $(CFLAGS) FileUtil.cpp

BamEvidence.o: BamEvidence.cpp BamEvidence.h
	$(CC) $(INCLUDES) $(CFLAGS) BamEvidence.cpp

EdgeEvidence.o: EdgeEvidence.cpp EdgeEvidence.h
	$(CC) $(INCLUDES) $(CFLAGS) EdgeEvidence.cpp

InsSizeEvidence.o: InsSizeEvidence.cpp InsSizeEvidence.h
	$(CC) $(INCLUDES) $(CFLAGS) InsSizeEvidence.cpp

MorphEvidence.o: MorphEvidence.cpp MorphEvidence.h
	$(CC) $(INCLUDES) $(CFLAGS) MorphEvidence.cpp

Fasta.o: Fasta.cpp Fasta.h
	$(CC) $(INCLUDES) $(CFLAGS) Fasta.cpp

clean:
	rm -rf *.o cnvthresher
