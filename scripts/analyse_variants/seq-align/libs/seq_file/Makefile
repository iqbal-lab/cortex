# To compile SAM/BAM support:
#     make HTSLIB=../htslib

ifdef HTSLIB
	ABSPATH=$(realpath $(HTSLIB))
	PARAMS := $(PARAMS) HTSLIB=$(ABSPATH)
endif

ifdef DEBUG
	PARAMS := $(PARAMS) DEBUG=$(DEBUG)
endif

#
# Compile dnacat
#
ifdef HTSLIB
	HTSARGS=-I $(HTSLIB) -D_USESAM=1 $(HTSLIB)/libhts.a -lm
endif

CFLAGS=-Wall -Wextra -std=c99 -pedantic -I.
LINKING=$(HTSARGS) -lpthread -lz

ifdef DEBUG
	OPT = -O0 --debug -g
else
	OPT = -O3
endif

all: bin/dnacat bin/dnademux benchmarks dev

bin/dnacat: tools/dna_cat.c seq_file.h stream_buffer.h
	mkdir -p bin
	$(CC) $(CFLAGS) $(OPT) -o $@ $< $(LINKING) -lm

bin/dnademux: tools/dna_demux.c seq_file.h stream_buffer.h
	mkdir -p bin
	$(CC) $(CFLAGS) $(OPT) -o $@ $< $(LINKING) -lm

benchmarks:
	cd benchmarks; make $(PARAMS)
dev:
	cd dev; make $(PARAMS)

clean:
	rm -rf bin
	cd benchmarks; make clean
	cd dev; make clean

.PHONY: all clean benchmarks dev
