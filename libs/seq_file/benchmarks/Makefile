ifdef HTSLIB
	HTSARGS=-I $(HTSLIB) -D_USESAM=1 $(HTSLIB)/libhts.a -lm
endif

CFLAGS=-Wall -Wextra -I ..
LINKING=$(HTSARGS) $(SAMLINK) -lpthread -lz

ifdef DEBUG
	OPT = -O0 --debug -g
else
	OPT = -O3
endif

all: seqtest ktest

seqtest: seq_test.c ../seq_file.h ../stream_buffer.h
	$(CC) $(OPT) $(CFLAGS) -o seqtest seq_test.c $(LINKING)
ktest: ktest.c
	$(CC) $(OPT) -Wall -Wextra -o ktest ktest.c -lz

clean:
	rm -rf seqtest ktest *.dSYM *.greg

.PHONY: all clean htslib
