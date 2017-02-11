CC ?= gcc

PLATFORM := $(shell uname)
COMPILER := $(shell ($(CC) -v 2>&1) | tr A-Z a-z )

ifdef DEBUG
	OPT = -O0 -DDEBUG=1 --debug -g -ggdb
else
	ifneq (,$(findstring clang,$(COMPILER)))
		OPT = -O3
	else
		OPT = -O4
		# TGTFLAGS = -fwhole-program
	endif
endif

CFLAGS = -Wall -Wextra -pedantic -std=c99 $(OPT)
OBJFLAGS = -fPIC
LIBFLAGS = -L. -lstrbuf -lz

all: libstrbuf.a strbuf_test

string_buffer.o: string_buffer.c string_buffer.h stream_buffer.h
	$(CC) $(CFLAGS) $(OBJFLAGS) -c string_buffer.c -o string_buffer.o

libstrbuf.a: string_buffer.o
	ar -csru libstrbuf.a string_buffer.o

strbuf_test: strbuf_test.c libstrbuf.a
	$(CC) $(CFLAGS) $(TGTFLAGS) strbuf_test.c -o strbuf_test $(LIBFLAGS)

test: strbuf_test
	./strbuf_test

clean:
	rm -rf string_buffer.o libstrbuf.a strbuf_test *.dSYM *.greg
	rm -rf tmp.strbuf.*.txt tmp.strbuf.*.txt.gz

.PHONY: all clean test

