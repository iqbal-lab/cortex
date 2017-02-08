/*
https://github.com/noporpoise/seq_file
Isaac Turner <turner.isaac@gmail.com>
Nov 2016, Public Domain
*/

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <getopt.h>

#include "seq_file.h"

const char *cmdstr = NULL;
FILE *fh1 = NULL, *fh2 = NULL;
gzFile gz1, gz2;
seq_format fmt = SEQ_FMT_UNKNOWN;

const char usage[] = "  Demultiplex input sequence\n"
"  -F,--fasta       print in FASTA format\n"
"  -Q,--fastq       print in FASTQ format\n"
"  -P,--plain       print in plain format\n"
"  -z,--gzip        gzip output\n";

const char shortopts[] = "hFQPz";
static struct option longopts[] =
{
  {"help",       no_argument,       NULL, 'h'},
  {"fasta",      no_argument,       NULL, 'F'},
  {"fastq",      no_argument,       NULL, 'Q'},
  {"plain",      no_argument,       NULL, 'P'},
  {"gzip",       no_argument,       NULL, 'z'},
  {NULL, 0, NULL, 0}
};


void print_usage(const char *err, ...)
__attribute__((noreturn))
__attribute__((format(printf, 1, 2)));

void print_usage(const char *err, ...)
{
  if(err != NULL) {
    fputc('\n', stderr);
    va_list argptr;
    fprintf(stderr, "Error: ");
    va_start(argptr, err);
    vfprintf(stderr, err, argptr);
    va_end(argptr);
    fputc('\n', stderr);
    fputc('\n', stderr);
  }

  fprintf(stderr, "Usage: %s [OPTIONS] <out1> <out2> [in]\n", cmdstr);
  fputs(usage, stderr);

  exit(EXIT_FAILURE);
}


static void print_read(const read_t *r, FILE *fh, gzFile gz)
{
  if(fh) {
    switch(fmt) {
      case SEQ_FMT_FASTA: seq_print_fasta(r, fh, 0); break;
      case SEQ_FMT_FASTQ: seq_print_fastq(r, fh, 0); break;
      case SEQ_FMT_PLAIN: fputs(r->seq.b, fh); fputc('\n', fh); break;
      default: fprintf(stderr, "Got value: %i\n", (int)fmt); exit(-1);
    }
  } else {
    switch(fmt) {
      case SEQ_FMT_FASTA: seq_gzprint_fasta(r, gz, 0); break;
      case SEQ_FMT_FASTQ: seq_gzprint_fastq(r, gz, 0); break;
      case SEQ_FMT_PLAIN: gzputs(gz, r->seq.b); gzputc(gz, '\n'); break;
      default: fprintf(stderr, "Got value: %i\n", (int)fmt); exit(-1);
    }
  }
}

int main(int argc, char **argv)
{
  cmdstr = argv[0];

  int fmt_set = 0;
  bool gzip_out = false;
  if(argc == 1) print_usage(NULL);

  // Arg parsing
  int c;
  opterr = 0; // silence getopt error messages

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': print_usage(NULL); break;
      case 'F': fmt_set++; fmt = SEQ_FMT_FASTA; break;
      case 'Q': fmt_set++; fmt = SEQ_FMT_FASTQ; break;
      case 'P': fmt_set++; fmt = SEQ_FMT_PLAIN; break;
      case 'z': gzip_out = true; break;
      case '?': /* BADCH getopt_long has already printed error */
        print_usage("Bad option: %s\n", argv[optind-1]);
      default: abort();
    }
  }

  size_t nargs = argc - optind;
  if(nargs < 2) print_usage("Need two output files");
  if(nargs > 3) print_usage("Can't have more than one input file");

  // open output files
  if(gzip_out) {
    gz1 = gzopen(argv[optind], "w");
    gz2 = gzopen(argv[optind+1], "w");
  } else {
    fh1 = fopen(argv[optind], "w");
    fh2 = fopen(argv[optind+1], "w");
  }

  seq_file_t *sf = seq_open(nargs == 3 ? argv[optind+2] : "-");

  read_t r1, r2;
  seq_read_alloc(&r1);
  seq_read_alloc(&r2);

  while(seq_read(sf, &r1)) {
    if(fmt == SEQ_FMT_UNKNOWN) { /* detect output format */
      if(seq_is_plain(sf)) fmt = SEQ_FMT_PLAIN;
      else if(seq_is_fasta(sf)) fmt = SEQ_FMT_FASTA;
      else fmt = SEQ_FMT_FASTQ;
    }
    print_read(&r1, fh1, gz1);
    if(!seq_read(sf, &r2)) { fprintf(stderr, "Odd number of reads\n"); }
    else print_read(&r2, fh2, gz2);
  }

  if(gzip_out) {
    gzclose(gz1);
    gzclose(gz2);
  } else {
    fclose(fh1);
    fclose(fh2);
  }

  seq_close(sf);

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);

  return EXIT_SUCCESS;
}
