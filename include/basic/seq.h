#ifndef STDSEQ_H_
#define STDSEQ_H_

#include <stdio.h>
#include <global.h>

typedef struct
{
  char *name;
  char *seq;  // sequence 
  char* qual; // qualities 
} Sequence;

//returns length of sequence read
int read_sequence_from_fasta(FILE*, Sequence * seq, int max_read_length);
int read_sequence_from_fastq(FILE*, Sequence * seq, int max_read_length);


void alloc_sequence(Sequence * seq, int max_read_length, int max_name_length);

void free_sequence(Sequence ** );

#endif /* STDSEQ_H_ */
