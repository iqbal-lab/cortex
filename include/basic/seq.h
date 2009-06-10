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
int read_full_entry_from_fasta(FILE*, Sequence * seq, int max_read_lengtn);
int read_sequence_from_fastq(FILE*, Sequence * seq, int max_read_length);

int read_sequence_from_fasta(FILE *fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry, int offset);

void alloc_sequence(Sequence * seq, int max_read_length, int max_name_length);

void free_sequence(Sequence ** );
void shift_last_kmer_to_start_of_sequence(Sequence * sequence, int length, short kmer_size);

#endif /* STDSEQ_H_ */
