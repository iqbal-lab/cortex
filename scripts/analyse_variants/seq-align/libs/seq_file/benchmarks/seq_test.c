#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include "seq_file.h"

#define DEFAULT_BUFSIZE (1<<20)

static struct option longopts[] =
{
// General options
  {"help",    no_argument, NULL, 'h'},
  {"no-buf",  no_argument, NULL, 'B'},
  {"no-zlib", no_argument, NULL, 'Z'},
  {NULL, 0, NULL, 0}
};

const char shortopts[] = "hBZ";

static void print_usage(const char *cmd)
{
  fprintf(stderr, "usage: %s [--no-buf|--no-zlib] <file>\n", cmd);
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
  bool use_buf = true, use_zlib = true;
  const char *path = NULL;
  int c;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': print_usage(argv[0]); break;
      case 'B': use_buf = false; break;
      case 'Z': use_zlib = false; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        print_usage(argv[0]);
      default: abort();
    }
  }

  if(optind >= argc) print_usage(argv[0]);
  path = argv[optind];

  seq_file_t *f = seq_open2(path, false, use_zlib, use_buf ? DEFAULT_BUFSIZE : 0);
  read_t r;
  seq_read_alloc(&r);
  if(f == NULL) { fprintf(stderr, "Cannot read: %s\n", path); exit(EXIT_FAILURE); }
  while(seq_read(f,&r) > 0)
    printf("%s\t[%lu,%lu,%lu]\n", r.name.b, r.name.end, r.seq.end, r.qual.end);
  seq_close(f);
  seq_read_dealloc(&r);
  return EXIT_SUCCESS;
}
