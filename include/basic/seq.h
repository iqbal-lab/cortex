#ifndef STDSEQ_H_
#define STDSEQ_H_

#include <stdio.h>

#define SEQ_MAX_NAME_LEN 255

#define INIT_SEQ(seq) (*seq).seq = NULL; (*seq).max=0;(*seq).name = NULL; (*seq).comment = NULL;

typedef struct
{
  char *name;
  char *comment;
  int length,max;
  char *seq; /* sequence */
  char* qual;
} Sequence;


Sequence * read_sequence_from_fasta(FILE*);
Sequence * read_sequence_from_fastq(FILE*);
void free_sequence(Sequence ** );

#endif /* STDSEQ_H_ */
