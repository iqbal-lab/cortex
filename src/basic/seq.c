
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <seq.h>

static int SEQ_BLOCK_SIZE = 512;

void seq_set_block_size(int size)
{
	SEQ_BLOCK_SIZE = size;
}

/* Read sequence from file "fp" in FASTA format.
   returns a pointer to seq (struct).
   it returns null if no sequence is left in file
*/

Sequence * read_sequence_from_fasta(FILE *fp)
{

// Note - If we do have a line bigger than 5000, fgets will read it in multiple consecutive goes.
//
  const int MAXLINE = 5000;
  char line[5000];
  int i;
  char name[256];


  if (fp == NULL)
  {
	  printf("Cannot pass a NULL pointer to read_sequence_from_fasta");
	  exit(1);
  }

  Sequence * seq = malloc(sizeof(Sequence));

  if (seq == NULL)
  {
	  printf("Out of memory allocating a Sequence object");
	  exit(1);
  }

  seq->seq = NULL;
  seq->max=0;
  seq->name = NULL;
  seq->comment = NULL;
  //seq->qual=NULL;  //uncomment this when you add qual's to struct seq.

  while (fgets(line, MAXLINE, fp) != NULL){
    //read name
	  // TODO fix it so it either works for > line longer than MAXLINE characters, or complains if that happens
	 	 // currently it won't work. fgets will get the excess of the line in the second loop, and it wont
	 	  // satisfy the if and so won't go into the name variable.

	  if (line[0] == '>')
	  {
      for(i = 1;i<MAXLINE;i++)
      {
    	  if (line[i] == '\n') //(line[i] == ' ' || line[i] == '\t' || line[i] == '\r')
    	  {
    		  break;
    	  }
    	  name[i-1] = line[i];
      }
      name[i-1] = '\0';

      seq->name = name;

      //read Sequence

      int length = 0;
      int max    = seq->max; //this is zero at the moment
      char c = 0;


      //Zam - either learn to understand this or change it so you do. TODO
      while (!feof(fp) && (c = fgetc(fp)) != '>') {
    	  if (isalpha(c) || c == '-' || c == '.') {
    		  if (length + 1 >= max) {
    			  max += SEQ_BLOCK_SIZE;
    			  seq->seq = (char*)realloc(seq->seq, sizeof(char) * max);
    		  }
    		  seq->seq[length++] = (char)c;
    	  }
      }
      if (c == '>') ungetc(c,fp);

      seq->seq[length] = '\0';
      seq->max         = max;
      seq->length      = length;

      return seq;
    }
    else{
      return NULL;
    }
  }

  return NULL;
}


Sequence * read_sequence_from_fastq(FILE *fp)
{

// Note - If we do have a line bigger than 5000, fgets will read it in multiple consecutive goes.
//
  const int MAXLINE = 5000;
  char line[5000];
  int i;
  char name[256];


  if (fp == NULL)
  {
	  printf("Cannot pass a NULL pointer to read_sequence_from_fastq");
	  exit(1);
  }

  Sequence * seq = malloc(sizeof(Sequence));

  if (seq == NULL)
  {
	  printf("Out of memory allocating a Sequence object");
	  exit(1);
  }

  seq->seq = NULL;
  seq->max=0;
  seq->name = NULL;
  seq->comment = NULL;
  seq->qual=NULL;  

  while (fgets(line, MAXLINE, fp) != NULL){
    //read name
	  // TODO fix it so it either works for > line longer than MAXLINE characters, or complains if that happens
	 	 // currently it won't work. fgets will get the excess of the line in the second loop, and it wont
	 	  // satisfy the if and so won't go into the name variable.

	  if (line[0] == '>')
	  {
      for(i = 1;i<MAXLINE;i++)
      {
    	  if (line[i] == '\n') //(line[i] == ' ' || line[i] == '\t' || line[i] == '\r')
    	  {
    		  break;
    	  }
    	  name[i-1] = line[i];
      }
      name[i-1] = '\0';

      seq->name = name;

      //read Sequence

      int length = 0;
      int max    = seq->max; //this is zero at the moment
      char c = 0;


      while (!feof(fp) && (c = fgetc(fp)) != '\n') {
    	  if (isalpha(c) || c == '-' || c == '.') {
    		  if (length + 1 >= max) {
    			  max += SEQ_BLOCK_SIZE;
    			  seq->seq = (char*)realloc(seq->seq, sizeof(char) * max);
    		  }
    		  seq->seq[length++] = (char)c;
    	  }
      }

      //ignore strand
      char* temp;
      fgets(fp,temp,2);

      //get quality
      int qual_length=0;
      int qual_max=0;

      while (!feof(fp) && (c = fgetc(fp)) != '\n') {
	if (length + 1 >= qual_max) {
	  qual_max += SEQ_BLOCK_SIZE;
	  seq->qual = (char*)realloc(seq->seq, sizeof(char) * qual_max);
	}
	seq->qual[qual_length++] = (char)c;
	
      }

      if ( (qual_max != max) || (qual_length != length))
	{
	  printf("Not reading fastq correctly");
	  exit(1);
	}

      seq->seq[length] = '\0';
      seq->qual[length]='\0';
      seq->max         = max;
      seq->length      = length;

      return seq;
    }
    else{
      return NULL;
    }
  }

  return NULL;
}




void free_sequence(Sequence ** sequence){
  free((*sequence)->seq);
  free((*sequence)->qual);
  free(*sequence);
  *sequence = NULL;
}

