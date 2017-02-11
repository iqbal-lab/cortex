/*
https://github.com/noporpoise/seq_file
Isaac Turner <turner.isaac@gmail.com>
Jan 2014, Public Domain
*/

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
#include <getopt.h>
#include <inttypes.h>
#include <ctype.h> // toupper() tolower() isprint()

#include <time.h>
#include <sys/time.h> // for seeding random
#include <unistd.h> // getpid()

#include "seq_file.h"

#define OPS_UPPERCASE       1 /* Convert to uppercase */
#define OPS_LOWERCASE       2 /* Convert to lowercase */
#define OPS_REVERSE         4 /* Reverse reads */
#define OPS_COMPLEMENT      8 /* Complement bases */
#define OPS_MASK_LC        16 /* Mask lower case bases */
#define OPS_NAME_ONLY      32 /* Print read name */
#define OPS_PRINT_LENGTH   64 /* Print read length */
#define OPS_KEY           128 /* Print lexically lower of seq and rev. complement */

const char usage[] = "  Read and manipulate dna sequence.\n"
#ifdef _USESAM
"  Compiled with SAM/BAM support.\n"
#endif
"\n"
"  -h,--help        show this help text\n"
"  -F,--fasta       print in FASTA format\n"
"  -Q,--fastq       print in FASTQ format\n"
"  -P,--plain       print in plain format\n"
"  -w,--wrap <n>    wrap lines by <n> characters [default: 0 (off)]\n"
"  -u,--uppercase   convert sequence to uppercase\n"
"  -l,--lowercase   convert sequence to lowercase\n"
"  -r,--revcmp      reverse complement sequence [i.e. -R and -C]\n"
"  -R,--reverse     reverse sequence\n"
"  -C,--complement  complement sequence\n"
"  -k,--key         give lexically lower of sequence and reverse complement\n"
"  -i,--interleave  interleave input files\n"
"  -m,--mask        mask lowercase bases\n"
"  -n,--rand <n>    print <n> random bases AFTER reading files\n"
"  -N,--names       print read names only\n"
"  -L,--lengths     print read names and lengths (implies -N)\n"
"  -s,--stat        probe and print file info, summarise read lengths\n"
"  -S,--fast-stat   probe and print file info only\n"
"  -M,--rename <f>  read names from <f>, one per line\n"
"\n"
"  Written by Isaac Turner <turner.isaac@gmail.com>\n";

static struct option longopts[] =
{
  {"help",       no_argument,       NULL, 'h'},
  {"fasta",      no_argument,       NULL, 'F'},
  {"fastq",      no_argument,       NULL, 'Q'},
  {"plain",      no_argument,       NULL, 'P'},
  {"wrap",       required_argument, NULL, 'w'},
  {"uppercase",  no_argument,       NULL, 'u'},
  {"lowercase",  no_argument,       NULL, 'l'},
  {"revcmp",     no_argument,       NULL, 'r'},
  {"reverse",    no_argument,       NULL, 'R'},
  {"complement", no_argument,       NULL, 'C'},
  {"interleave", no_argument,       NULL, 'i'},
  {"key",        no_argument,       NULL, 'k'},
  {"mask",       no_argument,       NULL, 'm'},
  {"rand",       required_argument, NULL, 'n'},
  {"names",      no_argument,       NULL, 'N'},
  {"lengths",    no_argument,       NULL, 'L'},
  {"stat",       no_argument,       NULL, 's'},
  {"fast-stat",  no_argument,       NULL, 'S'},
  {"rename",     required_argument, NULL, 'M'},
  {NULL, 0, NULL, 0}
};

const char shortopts[] = "hFQPw:ulrRCimn:NLsSM:";
const char *cmdstr;

const char bases[] = "ACGT";

#define die(fmt,...) do { \
  fprintf(stderr, "[%s:%i] Error: %s() "fmt"\n", __FILE__, __LINE__, __func__, __VA_ARGS__); \
  exit(EXIT_FAILURE); \
} while(0)

#define inpathstr(p) (strcmp(p,"-") == 0 ? "STDIN" : p)

char parse_entire_size(const char *str, size_t *result)
{
  char *strtol_last_char_ptr = NULL;
  if(*str < '0' || *str > '9') return 0;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);
  if(tmp > SIZE_MAX) return 0;
  if(strtol_last_char_ptr == NULL || *strtol_last_char_ptr != '\0') return 0;
  *result = (size_t)tmp;
  return 1;
}

size_t num_of_digits(size_t num)
{
  size_t digits = 1;
  while(1) {
    if(num < 10) return digits;
    if(num < 100) return digits+1;
    if(num < 1000) return digits+2;
    if(num < 10000) return digits+3;
    num /= 10000;
    digits += 4;
  }
  return digits;
}

// result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('18,446,744,073,709,551,615')+1 = 27
// returns pointer to result
char* ulong_to_str(unsigned long num, char* result)
{
  unsigned int digits = num_of_digits(num);
  unsigned int i, num_commas = (digits-1) / 3;
  char *p = result + digits + num_commas;
  *(p--) = '\0';

  for(i = 0; i < digits; i++, num /= 10) {
    if(i > 0 && i % 3 == 0) *(p--) = ',';
    *(p--) = '0' + (num % 10);
  }

  return result;
}

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

  fprintf(stderr, "Usage: %s [OPTIONS] <file1> [file2] ..\n", cmdstr);
  fputs(usage, stderr);

  exit(EXIT_FAILURE);
}

// 2 ops per byte h = strhash_fast_mix(h,x)
#define strhash_fast_mix(h,x) ((h) * 37 + (x))

static void seed_random()
{
  struct timeval now;
  gettimeofday(&now, NULL);

  uint32_t h;
  h = strhash_fast_mix(0, (uint32_t)now.tv_sec);
  h = strhash_fast_mix(h, (uint32_t)now.tv_usec);
  h = strhash_fast_mix(h, (uint32_t)getpid());
  srand(h);
}

// compare seq with the reverse complement of itself
int dna_rc_ncasecmp(const char *seq, size_t len)
{
  size_t i, j;
  for(i = 0, j = len-1; i < len; i++, j--) {
    int cmp = (int)tolower(seq[i]) - tolower(seq_char_complement(seq[j]));
    if(cmp) return cmp;
  }
  return 0;
}

// compare seq with the reverse of itself
int dna_r_ncasecmp(const char *seq, size_t len)
{
  size_t i, j;
  for(i = 0, j = len-1; i < len; i++, j--) {
    int cmp = (int)tolower(seq[i]) - tolower(seq[j]);
    if(cmp) return cmp;
  }
  return 0;
}

// compare seq with the complement of itself
int dna_c_ncasecmp(const char *seq, size_t len)
{
  size_t i;
  for(i = 0; i < len; i++) {
    int cmp = (int)tolower(seq[i]) - tolower(seq_char_complement(seq[i]));
    if(cmp) return cmp;
  }
  return 0;
}

static inline bool _print_rename_hdr(FILE *rename_fh, seq_buf_t *rnbuf,
                                     seq_format fmt)
{
  if((fmt == SEQ_FMT_FASTA || fmt == SEQ_FMT_FASTQ) && rename_fh)
  {
    rnbuf->end = 0; rnbuf->b[0] = '\0';
    if(!freadline(rename_fh, &rnbuf->b, &rnbuf->end, &rnbuf->size)) return false;
    cbuf_chomp(rnbuf->b, &rnbuf->end);

    // Strip off @ or > char at beginning
    if(rnbuf->end && (rnbuf->b[0] == '>' || rnbuf->b[0] == '@')) {
      memmove(rnbuf->b, rnbuf->b+1, --rnbuf->end);
      rnbuf->b[rnbuf->end] = '\0';
    }

    return true;
  }

  return false;
}

static void process_read(read_t *r, uint8_t ops)
{
  size_t i;

  if(ops & OPS_MASK_LC) {
    for(i = 0; i < r->seq.end; i++) {
      switch(r->seq.b[i]) {
        case 'A': case 'C': case 'G': case 'T': break;
        default: r->seq.b[i] = 'N';
      }
    }
  }

  if(ops & OPS_UPPERCASE)       seq_read_to_uppercase(r);
  else if(ops & OPS_LOWERCASE)  seq_read_to_lowercase(r);

  if((ops & OPS_REVERSE) && (ops & OPS_COMPLEMENT)) {
    if(!(ops & OPS_KEY) || dna_rc_ncasecmp(r->seq.b, r->seq.end) > 0)
      seq_read_reverse_complement(r);
  }
  else if((ops & OPS_REVERSE)) {
    if(!(ops & OPS_KEY) || dna_r_ncasecmp(r->seq.b, r->seq.end) > 0)
      seq_read_reverse(r);
  } else if((ops & OPS_COMPLEMENT)) {
    if(!(ops & OPS_KEY) || dna_c_ncasecmp(r->seq.b, r->seq.end) > 0)
      seq_read_complement(r);
  }
}

// Returns format used
static seq_format read_print(seq_file_t *sf, read_t *r,
                             seq_format fmt, uint8_t ops, size_t linewrap,
                             FILE *rename_fh, seq_buf_t *rnbuf)
{
  process_read(r, ops);

  if(fmt == SEQ_FMT_UNKNOWN) {
    // default to plain format is printing names only with no fmt specified
    if(ops & OPS_NAME_ONLY) fmt = SEQ_FMT_PLAIN;
    else if(seq_is_plain(sf)) fmt = SEQ_FMT_PLAIN;
    else if(seq_is_fasta(sf)) fmt = SEQ_FMT_FASTA;
    else fmt = SEQ_FMT_FASTQ;
  }

  /* Overwrite read name with rename buffer */
  if(_print_rename_hdr(rename_fh, rnbuf, fmt)) {
    // seq_buf_t tmp = r->name; r->name = *rnbuf; *rnbuf = tmp;
    cbuf_capacity(&r->name.b, &r->name.size, rnbuf->end);
    memcpy(r->name.b, rnbuf->b, rnbuf->end);
    r->name.b[r->name.end = rnbuf->end] = '\0';
  }

  if(ops & OPS_NAME_ONLY) {
    switch(fmt) {
      case SEQ_FMT_FASTA: fputc('>', stdout); break;
      case SEQ_FMT_FASTQ: fputc('@', stdout); break;
      case SEQ_FMT_PLAIN: break;
      default: die("Got value: %i\n", (int)fmt);
    }
    fputs(r->name.b, stdout);
    if(ops & OPS_PRINT_LENGTH) printf("\t%zu", r->seq.end);
    fputc('\n', stdout);
  }
  else {
    switch(fmt) {
      case SEQ_FMT_FASTA: seq_print_fasta(r, stdout, linewrap); break;
      case SEQ_FMT_FASTQ: seq_print_fastq(r, stdout, linewrap); break;
      case SEQ_FMT_PLAIN: fputs(r->seq.b, stdout); fputc('\n', stdout); break;
      default: die("Got value: %i\n", (int)fmt);
    }
  }

  return fmt;
}

static inline void _print_rnd_entries(const size_t *lens, size_t nentries,
                                      uint8_t fmt, size_t linewrap,
                                      FILE *rename_fh, seq_buf_t *rnbuf)
{
  size_t i, j, k, rnd = 0;

  for(i = 0; i < nentries; i++)
  {
    if(_print_rename_hdr(rename_fh, rnbuf, fmt))
      printf("%c%s\n", fmt == SEQ_FMT_FASTQ ? '@' : '>', rnbuf->b);
    else if(fmt == SEQ_FMT_FASTA) printf(">rand%zu\n", i);
    else if(fmt == SEQ_FMT_FASTQ) printf("@rand%zu\n", i);

    for(j = k = 0; j < lens[i]; j++, k++) {
      if(linewrap && k == linewrap) { k = 0; fputc('\n', stdout); }
      // use 2 bits per iteration, 32 bits in rand(), update every 16 iterations
      if((j & 15) == 0) rnd = (size_t)rand();
      fputc(bases[rnd&3], stdout);
      rnd >>= 2;
    }
    if(fmt == SEQ_FMT_FASTQ) { /* quality scores */
      fputs("\n+\n", stdout);
      for(j = k = 0; j < lens[i]; j++, k++) {
        if(linewrap && k == linewrap) { k = 0; fputc('\n', stdout); }
        fputc(33+rand()%41, stdout); // 33..73
      }
    }
    fputc('\n', stdout);
  }
}

// @fast if true skip reading over all reads
static void file_stat(seq_file_t *sf, read_t *r, uint8_t ops, bool fast)
{
  if(ops && fast)
    print_usage("-S,--fast-stat is not compatible with -l,-u,-r,-R,-C,-m");

  if(ops & ~(OPS_UPPERCASE | OPS_LOWERCASE))
    print_usage("-s,--stat and -S,--fast-stat are not compatible with -r,-R,-C,-m");

  printf("File: %s\n", inpathstr(sf->path));

  int minq = -1, maxq = -1, s, fmti;

  fmti = seq_guess_fastq_format(sf, &minq, &maxq);
  s = seq_read(sf,r);

  if(s < 0) die("Error reading file: %s\n", inpathstr(sf->path));
  if(s == 0) die("Cannot get any reads from file: %s\n", inpathstr(sf->path));

  if(seq_is_sam(sf)) printf("  Format: SAM\n");
  if(seq_is_bam(sf)) printf("  Format: BAM\n");
  if(seq_is_fasta(sf)) printf("  Format: FASTA\n");
  if(seq_is_fastq(sf)) printf("  Format: FASTQ\n");
  if(seq_is_plain(sf)) printf("  Format: plain\n");

  if(seq_use_gzip(sf)) printf("  Read with zlib\n");

  char print_qstat = (seq_is_fastq(sf) || seq_is_sam(sf) || seq_is_bam(sf));

  if(print_qstat)
  {
    if(fmti == -1) printf("  Couldn't get any quality scores\n");
    else {
      printf("  Format QScores: %s, offset: %i, min: %i, max: %i, scores: [%i,%i]\n",
             FASTQ_FORMATS[fmti], FASTQ_OFFSET[fmti], FASTQ_MIN[fmti], FASTQ_MAX[fmti],
             FASTQ_MIN[fmti]-FASTQ_OFFSET[fmti], FASTQ_MAX[fmti]-FASTQ_OFFSET[fmti]);
      printf("  QScore range in first 500bp: [%i,%i]\n", minq, maxq);
    }
  }

  if(!fast)
  {
    // We've already read one read
    size_t i, char_count[256] = {0};
    size_t total_len = 0, nreads = 0;
    size_t min_rlen = SIZE_MAX, max_rlen = 0;

    do {
      process_read(r, ops);
      total_len += r->seq.end;
      max_rlen = r->seq.end > max_rlen ? r->seq.end : max_rlen;
      min_rlen = r->seq.end < min_rlen ? r->seq.end : min_rlen;
      nreads++;

      for(i = 0; i < r->seq.end; i++)
        char_count[(uint8_t)r->seq.b[i]]++;
    }
    while((s = seq_read(sf,r)) > 0);

    size_t mean_rlen = (size_t)(((double)total_len / nreads) + 0.5);

    char nbasesstr[50], nreadsstr[50];
    char minrlenstr[50], maxrlenstr[50], meanrlenstr[50];
    ulong_to_str(total_len, nbasesstr);
    ulong_to_str(nreads, nreadsstr);
    ulong_to_str(min_rlen, minrlenstr);
    ulong_to_str(max_rlen, maxrlenstr);
    ulong_to_str(mean_rlen, meanrlenstr);

    printf("  Total seq (bp):     %s\n", nbasesstr);
    printf("  Number of reads:    %s\n", nreadsstr);
    printf("  Shortest read (bp): %s\n", minrlenstr);
    printf("  Longest read  (bp): %s\n", maxrlenstr);
    printf("  Mean length   (bp): %s\n", meanrlenstr);

    printf("  Char Counts:\n");
    for(i = 0; i < 256; i++) {
      if(char_count[i]) {
        if(isprint((char)i)) printf("      %c: %zu\n", (char)i, char_count[i]);
        else                 printf(" (%3zu): %zu\n",        i, char_count[i]);
      }
    }

    if(s < 0) die("Error reading file: %s\n", inpathstr(sf->path));
  }

  fputc('\n', stdout);
}

static void vector_push(size_t **ptr, size_t *len, size_t *cap, size_t x)
{
  if(!*ptr || *len >= *cap) {
    if(!*cap) { *cap = 16; }
    else { while(*len >= *cap) { *cap *= 2; }}
    *ptr = realloc(*ptr, *cap * sizeof(size_t));
  }
  (*ptr)[(*len)++] = x;
}

int main(int argc, char **argv)
{
  cmdstr = argv[0];

  bool interleave = false, stat = false, fast_stat = false;
  uint8_t ops = 0, fmt_set = 0;
  seq_format fmt = SEQ_FMT_UNKNOWN;
  size_t i, linewrap = 0;
  char *rename_path = NULL;

  size_t *nrand = NULL, nrand_len = 0, nrand_cap = 0, tmprnd = 0;

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
      case 'w':
        if(!parse_entire_size(optarg, &linewrap))
          print_usage("Bad -w argument: %s\n", optarg);
        break;
      case 'u': ops |= OPS_UPPERCASE;    break;
      case 'l': ops |= OPS_LOWERCASE;    break;
      case 'r': ops |= OPS_REVERSE | OPS_COMPLEMENT; break;
      case 'R': ops |= OPS_REVERSE;      break;
      case 'C': ops |= OPS_COMPLEMENT;   break;
      case 'm': ops |= OPS_MASK_LC;      break;
      case 'N': ops |= OPS_NAME_ONLY;    break;
      case 'L': ops |= OPS_NAME_ONLY | OPS_PRINT_LENGTH; break;
      case 'k': ops |= OPS_KEY;          break;
      case 'n':
        if(!parse_entire_size(optarg, &tmprnd))
          print_usage("Bad -n argument: %s\n", optarg);
        vector_push(&nrand, &nrand_len, &nrand_cap, tmprnd);
        break;
      case 'i': interleave  = true;   break;
      case 's': stat        = true;   break;
      case 'S': fast_stat   = true;   break;
      case 'M': rename_path = optarg; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        print_usage("Bad option: %s\n", argv[optind-1]);
      default: abort();
    }
  }

  if(fmt_set > 1)
    print_usage("Please specify only one output format (-f,-q,-p)\n");

  size_t num_inputs = argc - optind;
  char **input_paths = argv + optind;

  if(!nrand_len && !num_inputs)
    print_usage("Please specify at least one input file\n");

  // Default to plain format for random output
  if(nrand_len && !num_inputs && fmt == SEQ_FMT_UNKNOWN) fmt = SEQ_FMT_PLAIN;

  if(linewrap && (fmt == SEQ_FMT_PLAIN))
    print_usage("Bad idea to use linewrap with plain output (specify -f or -q)");

  if((ops & OPS_UPPERCASE) && (ops & OPS_LOWERCASE))
    print_usage("Cannot use both -u,--uppercase and -l,--lowercase");

  if((ops & OPS_NAME_ONLY) &&
     ((ops &~(OPS_NAME_ONLY|OPS_PRINT_LENGTH)) || nrand_len ||
      (fmt != SEQ_FMT_PLAIN && fmt != SEQ_FMT_UNKNOWN))) {
    print_usage("Cannot use -N,--names or -L,--lengths with other options");
  }

  if((ops & OPS_KEY) && !(ops & (OPS_REVERSE | OPS_COMPLEMENT)))
    print_usage("Must use --key with one of --reverse, --complment, --revcmp");

  if(stat && (interleave || linewrap || fmt || nrand_len || ops))
    print_usage("-s,--stat is not compatible with other options");

  if(stat && fast_stat)
    print_usage("Cannot use -s,--stat and -S--fast-stat together");

  FILE *rename_fh = NULL;
  seq_buf_t rename_buf;

  if(rename_path) {
    if(strcmp(rename_path,"-") == 0) rename_fh = stdin;
    else if((rename_fh = fopen(rename_path, "r")) == NULL)
      die("Cannot open --rename file: %s", inpathstr(rename_path));
    rename_buf.size = 1024;
    rename_buf.b = malloc(rename_buf.size);
    if(!rename_buf.b) die("Out of memory%c", '!');
  }

  if(nrand_len) seed_random();

  read_t r;
  seq_read_alloc(&r);
  int s;

  seq_file_t *inputs[num_inputs];

  for(i = 0; i < num_inputs; i++) {
    if((inputs[i] = seq_open(input_paths[i])) == NULL)
      print_usage("Couldn't read file: %s\n", inpathstr(input_paths[i]));
  }

  if(stat || fast_stat) {
    for(i = 0; i < num_inputs; i++)
      file_stat(inputs[i], &r, ops, fast_stat);
  }
  else if(interleave) {
    // read one entry from each file
    size_t waiting_files = num_inputs;
    while(waiting_files) {
      for(i = 0; i < num_inputs; i++) {
        if(inputs[i] != NULL) {
          s = seq_read(inputs[i],&r);
          if(s < 0) die("Error reading from: %s\n", inputs[i]->path);
          else if(s > 0) {
            fmt = read_print(inputs[i], &r, fmt, ops, linewrap,
                             rename_fh, &rename_buf);
          } else {
            seq_close(inputs[i]); inputs[i] = NULL; waiting_files--;
          }
        }
      }
    }
  }
  else {
    for(i = 0; i < num_inputs; i++) {
      while((s = seq_read(inputs[i],&r)) > 0) {
        fmt = read_print(inputs[i], &r, fmt, ops, linewrap,
                         rename_fh, &rename_buf);
      }
      if(s < 0) die("Error reading from: %s\n", inputs[i]->path);
      seq_close(inputs[i]);
    }
  }

  seq_read_dealloc(&r);

  // Print random entries
  _print_rnd_entries(nrand, nrand_len, fmt, linewrap, rename_fh, &rename_buf);
  free(nrand);

  if(rename_fh) {
    fclose(rename_fh);
    free(rename_buf.b);
  }

  return EXIT_SUCCESS;
}
