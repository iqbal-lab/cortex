/*
  binary_kmer.h - basic manipulation of binary kmers
  
*/

#ifndef BINARY_KMER_H_
#define BINARY_KMER_H_
#include <global.h>

// *** Basic constants to allow compile-time configuration of what kmer-sizes we support ****
#define ALLOW_KMERS_UP_TO_31 1
#define ALLOW_KMERS_UP_TO_63 2
#define ALLOW_KMERS_UP_TO_95 3
#define ALLOW_KMERS_UP_TO_127 4
#define ALLOW_KMERS_UP_TO_160 5
#define ALLOW_KMERS_UP_TO_192 6
#define ALLOW_KMERS_UP_TO_223 7
#define ALLOW_KMERS_UP_TO_255 8

// change this if you want to support different kmer-ranges
#define NUMBER_OF_BITFIELDS_IN_BINARY_KMER  ALLOW_KMERS_UP_TO_31



typedef enum
 {
    Adenine   = 0,
    Cytosine  = 1,
    Guanine   = 2,
    Thymine   = 3,
    Undefined = 4,
  } Nucleotide ;


//platform specific
typedef unsigned long long bitfield_of_64bits;

typedef bitfield_of_64bits BinaryKmer[NUMBER_OF_BITFIELDS_IN_BINARY_KMER];  //think of this as the number of long longs we encode the kmer in



typedef struct
{
  int nkmers;
  BinaryKmer * kmer;
} KmerSlidingWindow;


//a set of KmerArrays
typedef struct
{
  int nwindows;
  int max_nwindows;
  KmerSlidingWindow * window; 
} KmerSlidingWindowSet;


// basic BinaryKmer operations
void       binary_kmer_initialise_to_zero(BinaryKmer* bkmer);
void       binary_kmer_assignment_operator(BinaryKmer left, BinaryKmer right);
boolean    binary_kmer_comparison_operator(const BinaryKmer const left, const BinaryKmer const right);
boolean    binary_kmer_less_than(const BinaryKmer const left, const BinaryKmer const right,short kmer_size);
void       binary_kmer_right_shift(BinaryKmer* kmer, int num_bits_to_shift);
void       binary_kmer_left_shift(BinaryKmer* kmer, int num_bits_to_shift, short kmer_size);
void       binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(BinaryKmer* bkmer, Nucleotide n, short kmer_size);



Nucleotide char_to_binary_nucleotide(char c);

char binary_nucleotide_to_char(Nucleotide n);

char * nucleotides_to_string(Nucleotide * nucleotides, int length, char * string);

Nucleotide reverse_binary_nucleotide(Nucleotide n);

//get overlapping kmers from sequence
int get_sliding_windows_from_sequence(char * sequence,char * qualities, int length, char quality_cutoff, short kmer_size, KmerSlidingWindowSet * windows, int max_windows, int max_kmers); 

//use preallocated sliding window, and get all the kmers from the passed-in sequence. Any kmer that would have contained an N is returned as NULL
int get_single_kmer_sliding_window_from_sequence(char * seq, int length, short kmer_size, KmerSlidingWindow* kmer_window);

char * binary_kmer_to_seq(BinaryKmer* kmer, short kmer_size, char * seq);

BinaryKmer* seq_to_binary_kmer(char * seq, short kmer_size, BinaryKmer* prealloced_kmer);

char * seq_reverse_complement(char * in, int length, char * out);

BinaryKmer* binary_kmer_reverse_complement(BinaryKmer* kmer, short kmer_size, BinaryKmer* prealloc_reverse_kmer);

Nucleotide binary_kmer_get_last_nucleotide(BinaryKmer* kmer);

Nucleotide binary_kmer_get_first_nucleotide(BinaryKmer* kmer, short kmer_size);

//BinaryKmer binary_kmer_add_nucleotide_shift(BinaryKmer kmer,Nucleotide nucleotide, short kmer_size);

char reverse_char_nucleotide(char c);

void binary_kmer_alloc_kmers_set(KmerSlidingWindowSet * windows, int max_windows, int max_kmers);
void binary_kmer_free_kmers(KmerSlidingWindow * *);
void binary_kmer_free_kmers_set(KmerSlidingWindowSet * *);

void nucleotide_iterator(void (*f)(Nucleotide));

#endif /* BINARY_KMER_H_ */
