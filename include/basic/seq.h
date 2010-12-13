#ifndef STDSEQ_H_
#define STDSEQ_H_

#include <stdio.h>
#include <global.h>

typedef struct
{
  char *name;
  int  start,end;
  char *seq;  // sequence 
  char* qual; // qualities 
} Sequence;

//returns length of sequence read

//note that fastq format dosn't support "partial" reading of a long entry -- ask Zam or Mario
//that means we can only read full fastq entries
int read_sequence_from_fastq(FILE * fp, Sequence * seq, int max_read_length, int fastq_ascii_offset);

//this routine can read long sequences (eg full chromosomes) , this is implemented by reading the sequence in chunks
int read_sequence_from_fasta(FILE * fp, Sequence * seq, int max_chunk_length, boolean new_entry, boolean * full_entry, int offset);

void alloc_sequence(Sequence * seq, int max_read_length, int max_name_length);

void free_sequence(Sequence ** );

void shift_last_kmer_to_start_of_sequence(Sequence * sequence, int length, short kmer_size);

#endif /* STDSEQ_H_ */
