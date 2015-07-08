seq_file
========

C Library for reading multiple bioinformatics sequence file formats  
https://github.com/noporpoise/seq_file  
Isaac Turner <turner.isaac@gmail.com>  
21 July 2014, license: Public Domain

[![Build Status](https://travis-ci.org/noporpoise/seq_file.png?branch=master)](https://travis-ci.org/noporpoise/seq_file)

About
=====

The aim is to provide a C library that allows programs to transparently read
sequence data without having to worry about the file format.
Pass a file to `seq_open(path)` and then read sequences using `seq_read(...)`.
File format is automatically detected.

Currently supports:
* SAM & BAM (if compiled with htslib with `make HTSLIB=`)
* FASTA (& gzipped fasta)
* FASTQ (& gzipped fastq)
* 'plain' format (one sequence per line [.txt]) (& gzipped plain [.txt.gz])

`seq_open_fh(...)` allows you to read through pipes and the command line.
`seq_open2(...)` gives more options about how you'd like to read your input.

Build
=====

seq_file requires [htslib](https://github.com/samtools/htslib) for reading SAM/BAM
files.  To compile run with htslib:

    make HTSLIB=PATH/TO/htslib

or without (you won't be able to read sam/bam files)

    make

Tools
=====

Also included is a tool for reading, converting amd manipulating sequence files: `dnacat`

Generate a random 20bp sequence and a 5bp sequence, print in fasta format

    ./bin/dnacat -n 20 -n 5 -F

Output:

    >rand0
    CTGATAAATCTTTTACGTCC
    >rand1
    GTAGC

Read three files, print as FASTQ:

    ./bin/dnacat -Q input.sam input.fq input.fa

Get some stats on an input file:

    ./bin/dnacat -s input.sam

Output:

    File: test/reads.fq
    Format: FASTQ (read with zlib)
    Format QScores: Sanger (Phred+33), offset: 33, min: 33, max: 73, scores: [0,40]
    QScore range in first 500bp: [35,65]
    Total seq (bp):     1,275
    Number of reads:    25
    Shortest read (bp): 51
    Longest read  (bp): 51
    Mean length   (bp): 51

Interleave two files, print in FASTQ:

    ./bin/dnacat -Q -i reads.1.fq.gz reads.2.fq.gz > reads.interleaved.fq.gz

Read a file, print in FASTA, wrap lines to be no longer than 80bp:

    ./bin/dnacat -F -w 80 input.fq > input.wrapped.fa

Other useful one liners:
* Reverse complement a sequence: `./bin/dnacat -r - <<< AACGA` Output: `TCGTT`
* Reverse a sequence: `./bin/dnacat -R - <<< AACGA` Output: `AGCAA`
* Complement a sequence: `./bin/dnacat -C - <<< AACGA` Output: `TTGCT`
* Convert to uppercase: `./bin/dnacat -u in.fa`
* Convert to lowercase: `./bin/dnacat -l in.fa`
* Convert lowercase + non-ACGT bases to 'N': `./bin/dnacat -m in.fa`

Note: in bash `cmd <<< AACGA` is the same as 'echo AACGA | cmd'

Example Code
============

Example code to read a file and print as a FASTA file using seq_file.

    #include "seq_file.h"

    int main(int argc, char **argv)
    {
      if(argc != 2) exit(EXIT_FAILURE);

      seq_file_t *file = seq_open(argv[1]);

      if(file == NULL)
        exit(EXIT_FAILURE);

      read_t *read = seq_read_new();

      while(seq_read(file, read) > 0)
      {
        printf("@%s\n", read->name.b);
        printf("%s\n", read->seq.b);
        printf("+\n);
        printf("%s\n", read->qual.b);
      }

      seq_close(file);
      seq_read_free(read);
    }

If the code above is pasted into a file test.c and the files `seq_file.h` and
`stream_buffer.h` are copied into the same directory, the program should compile
with the following command:

    gcc -Wall -o test test.c -lz

Or to use htslib:

    HTSLIB=PATH/TO/HTSLIB
    gcc -Wall -o test test.c -D_USESAM=1 -I$(HTSLIB)/htslib -L$(HTSLIB)/htslib -lhts -lpthread -lz

You'll need to compile htslib first.

Functions
=========

Opening/closing
---------------

    read_t* seq_read_new()

Returns a pointer to a new `read_t` struct

    void seq_read_free(read_t* r)

Free a `read_t` struct

    read_t* seq_read_alloc(read_t *r)

Allocate an object on the stack

    void seq_read_dealloc(read_t* r)

Free an object on the stack

    seq_file_t* seq_open(const char *path)

Open sequence file pointed to by path.

    seq_file_t* seq_open2(const char *path, char sam_bam, char use_gzip, size_t buffer_size)

Parameters:
* `sam_bam`: if 1 opens as sam, if 2 opens as bam, 0 treats as other
* `use_zlib`: if 0 opens with FILE, otherwise uses gzFile
* `buffer_size`: size of buffer to read into.

If `buffer_size == 0` no buffer is used, with the exception of `use_zlib=1` with newer versions of zlib where all input is buffered by the zlib library

    seq_file_t* seq_open_fh(FILE *fh, char buffered)

Use a file handle that has already been opened.  Pass `stdin` to read from cmdline,
pipes etc. If buffered is passed, attempts to read with a gzip stream and a buffer.
Note: seq_open_fh can only read sam/bam files from stdin at the moment.
`seq_close()` will close `fh` so you shouldn't call `fclose(fh)`. Returns `NULL` on error, in which case `fh` will not have been closed.

    void seq_close(seq_file_t *sf)

Close a `seq_file_t`

Reading
-------

    int seq_read(seq_file_t *sf, read_t *r)

Read a read from the file into `r`.
Returns 1 on success, 0 on eof, -1 if partially read / syntax error

    void seq_read_reverse_complement(read_t *r)

Reverse complement a read. If the read has quality scores, they are also reversed. 

File status
-----------

These functions return 1 (if true), otherwise 0

    seq_is_bam(seq_file_t *sf)
    seq_is_sam(seq_file_t *sf)
    seq_use_gzip(seq_file_t *sf) // is this seq_file reading through zlib?

The following require a read to have been read successfully using seq_read:

    seq_is_fasta(seq_file_t *sf)
    seq_is_fastq(seq_file_t *sf)
    seq_is_plain(seq_file_t *sf)

Probe a file
------------

    seqtype_t seq_guess_filetype_from_extension(const char *path)

Returns one of: {IS_SAM, IS_BAM, IS_FASTQ, IS_FASTA, IS_PLAIN, IS_UNKNOWN}

    int seq_guess_fastq_format(const char *path, int max_read_bases,
                               int *min_qual, int *max_qual)

Returns and number between 0 and 5 inclusive on success. These refer to:

0) Sanger / Illumina 1.9+ (Phred+33) -- this is the most common
1) Sanger (Phred+33)
2) Solexa (Solexa+64)
3) Illumina 1.5+ (Phred+64)
4) Illumina 1.3+ (Phred+64)
5) Illumina 1.8+ (Phred+33)

Returns -1 if it cannot read the file or the file has no quality scores (e.g
it's in FASTA format).

Example usage:

    int fmt = seq_guess_fastq_format("in.fq");
    printf("in.fq format: %s, offset: %i, min: %i, max: %i, range: [%i,%i]\n",
           FASTQ_FORMATS[fmt], FASTQ_MIN[fmt], FASTQ_MAX[fmt], FASTQ_OFFSET[fmt],
           FASTQ_MIN[fmt]-FASTQ_OFFSET[fmt], FASTQ_MAX[fmt]-FASTQ_OFFSET[fmt]);

Outputs:

    in.fq format: Sanger / Illumina 1.9+ (Phred+33), offset: 33, min: 126, max: 33, range: [0,93]
    in.fq format: Sanger (Phred+33), offset: 33, min: 73, max: 33, range: [0,40]
    in.fq format: Solexa (Solexa+64), offset: 59, min: 104, max: 64, range: [-5,40]
    in.fq format: Illumina 1.3+ (Phred+64), offset: 64, min: 104, max: 64, range: [0,40]
    in.fq format: Illumina 1.5+ (Phred+64), offset: 67, min: 104, max: 64, range: [3,40]
    in.fq format: Illumina 1.8+ (Phred+33), offset: 33, min: 74, max: 33, range: [0,41]

To get highest and lowest quality characters in the first `num` bases:

    int seq_get_qual_limits(const char *path, int num, int *minptr, int *maxptr)

`minptr` and `maxptr` are set to the min and max values of the first `num` base
quality scores.  Returns 0 if no qual scores in the first `num` bases,
1 on success, -1 if read error.

Writing
-------

    int seq_print_fasta(const read_t *r, FILE *fh, int linewrap)
    int seq_print_fastq(const read_t *r, FILE *fh, int linewrap)
    int seq_gzprint_fasta(const read_t *r, gzFile gz, int linewrap)
    int seq_gzprint_fastq(const read_t *r, gzFile gz, int linewrap)

Write a read in FASTA or FASTQ format.  If using FASTQ and the quality score
and sequence have different lengths, the quality score is shortened or padded
with '.' to make it the same length as the sequence. Return -1 on error, 0 on
success.

Useful functions
----------------

    static inline size_t seq_read_truncate_name(read_t *r)

Formally, FASTA/Q entry names stop at the first space character. 
This function truncates the read name and returns new length

    static inline void seq_read_to_uppercase(read_t *r)
    static inline void seq_read_to_lowercase(read_t *r)

Convert read sequence to upper or lower case.  Does not affect read name or
quality scores.

    seq_read_looks_valid_dna(read_t *r)
    seq_read_looks_valid_rna(read_t *r)
    seq_read_looks_valid_protein(read_t *r)

Check that a read looks valid.  Returns 1 if it appears valid, 0 otherwise.
A read is *invalid* if it has:
* an invalid sequence character
* quality scores but of a different length than the sequence
* a quality score that is `<33` or `>105`

Valid sequence characters are (upper and lower case are valid):
* DNA: ACGTN
* RNA: ACGUN
* protein: ACDEFGHIKLMNOPQRSTUVWY

License
=======

Public Domain. Use as you wish. No warranty. There may be bugs.
