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

all: seqcmdline

seqcmdline: seq_cmdline.c ../seq_file.h ../stream_buffer.h
	$(CC) $(OPT) $(CFLAGS) -o seqcmdline seq_cmdline.c $(LINKING)

clean:
	rm -rf seqcmdline *.dSYM *.greg

.PHONY: all clean
