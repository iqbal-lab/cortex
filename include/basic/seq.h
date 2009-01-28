#ifndef STDSEQ_H_
#define STDSEQ_H_

#include <stdio.h>


typedef struct
{
  char *name;
  char *comment;
  char *seq;  // sequence 
  char* qual; // qualities 
} Sequence;


int read_sequence_from_fasta(FILE*, Sequence * seq, int max_read_length);
int read_sequence_from_fastq(FILE*, Sequence * seq, int max_read_length);

void free_sequence(Sequence ** );

#endif /* STDSEQ_H_ */
