/*
  binary_kmer.h - basic manipulation of binary kmers
  
*/

#ifndef BINARY_KMER_H_
#define BINARY_KMER_H_

typedef enum
 {
    Adenine  = 0,
    Cytosine = 1,
    Guanine  = 2,
    Thymine  = 3,
    Undefined = 4
  } Nucleotide ;

typedef long long BinaryKmer;


// KmerArray is an array of kmers. Usually this will be from a sliding window run across a sequence.

typedef struct{
	int nkmers;
	BinaryKmer * bin_kmers;
} KmerArray;


Nucleotide char_to_binary_nucleotide(char c);

char binary_nucleotide_to_char(Nucleotide n);

//get overlapping kmers from sequence
KmerArray * get_binary_kmers_from_sequence(char * sequence,int sequence_length, short kmer_size); 

char * binary_kmer_to_seq(BinaryKmer kmer, short kmer_size);

BinaryKmer seq_to_binary_kmer(char * seq, short kmer_size);

BinaryKmer binary_kmer_reverse_complement(BinaryKmer kmer, short kmer_size);

Nucleotide binary_kmer_get_last_nucleotide(BinaryKmer kmer);

Nucleotide binary_kmer_get_first_nucleotide(BinaryKmer kmer, short kmer_size);

BinaryKmer binary_kmer_add_nucleotide_shift(BinaryKmer kmer,Nucleotide nucleotide, short kmer_size);

char reverse_char_nucleotide(char c);

void binary_kmer_free_kmers(KmerArray * *);

#endif /* BINARY_KMER_H_ */
