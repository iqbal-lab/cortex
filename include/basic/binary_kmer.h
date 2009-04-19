/*
  binary_kmer.h - basic manipulation of binary kmers
  
*/

#ifndef BINARY_KMER_H_
#define BINARY_KMER_H_

typedef enum
 {
    Adenine   = 0,
    Cytosine  = 1,
    Guanine   = 2,
    Thymine   = 3,
    Undefined = 4
  } Nucleotide ;

typedef long long BinaryKmer;


// KmerArray is an array of kmers. Usually this will be from a sliding window run across a sequence.
typedef struct{
	int nkmers;
	BinaryKmer * kmer;
} KmerSlidingWindow;


//a set of KmerArrays
typedef struct{
	int nwindows;
	int max_nwindows;
	KmerSlidingWindow * window; 
} KmerSlidingWindowSet;


Nucleotide char_to_binary_nucleotide(char c);

char binary_nucleotide_to_char(Nucleotide n);

char * nucleotides_to_string(Nucleotide * nucleotides, int length, char * string);

Nucleotide reverse_binary_nucleotide(Nucleotide n);

//get overlapping kmers from sequence
int get_sliding_windows_from_sequence(char * sequence,char * qualities, int length, char quality_cutoff, short kmer_size, KmerSlidingWindowSet * windows, int max_windows, int max_kmers); 

char * binary_kmer_to_seq(BinaryKmer kmer, short kmer_size, char * seq);

BinaryKmer seq_to_binary_kmer(char * seq, short kmer_size);

BinaryKmer binary_kmer_reverse_complement(BinaryKmer kmer, short kmer_size);

Nucleotide binary_kmer_get_last_nucleotide(BinaryKmer kmer);

Nucleotide binary_kmer_get_first_nucleotide(BinaryKmer kmer, short kmer_size);

BinaryKmer binary_kmer_add_nucleotide_shift(BinaryKmer kmer,Nucleotide nucleotide, short kmer_size);

char reverse_char_nucleotide(char c);

void binary_kmer_alloc_kmers_set(KmerSlidingWindowSet * windows, int max_windows, int max_kmers);
void binary_kmer_free_kmers(KmerSlidingWindow * *);
void binary_kmer_free_kmers_set(KmerSlidingWindowSet * *);

void nucleotide_iterator(void (*f)(Nucleotide));

#endif /* BINARY_KMER_H_ */
