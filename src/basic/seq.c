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


      while (!feof(fp) && (c = fgetc(fp)) != '>') {
    	  if (isalpha(c) || c == '-' || c == '.') {
    		  if (length + 1 >= max) {
    			  max += SEQ_BLOCK_SIZE;
			  seq->seq = (char*)realloc(seq->seq, sizeof(char) * max);
			  if (seq->seq == NULL)
			    {
			      printf("Out of mem doing realloc of seq");
			      exit(1);
			    }
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
      free(seq);
      return NULL;
    }
  }
  
  free(seq);
  return NULL;
}


//if going to be reading lots, alloc your own Sequence* just once,
// and call this repeatedly, reusing that memory.
//returns 0 in case of error. sequence is returned in the seq* you pass in as argument, which you must prealloc, including seq->seq.
int  read_sequence_from_fasta_efficient(FILE *fp, Sequence* seq, int* did_a_realloc_occur)
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


  if (seq == NULL)
  {
	  printf("Dont give me a null pointer for seq - alloc memory yourself and give me that");
	  exit(1);
  }
  if (seq->seq == NULL)
  {
	  printf("Dont give me a null pointer for seq->seq - alloc memory yourself and give me that");
	  exit(1);
  }



  //seq->seq = NULL;
  //  seq->max=0;
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
	int max    = seq->max; //This starts off as whatever size you alloced the seq->seq that you passed in
	char c = 0;
	
	
	while (!feof(fp) && (c = fgetc(fp)) != '>') {
    	  if (isalpha(c)) 
	    {
	      if (length + 1 >= max) 
		{
		  max += SEQ_BLOCK_SIZE;
		  seq->seq = (char*)realloc(seq->seq, sizeof(char) * max);
		  printf("Realloc of seq->seq hppening. I don't expect this to ever happen");
		  *did_a_realloc_occur=1;
		  if (seq->seq == NULL)
		    {
		      printf("Out of mem doing realloc of seq");
		      exit(1);
		    }
		}
	      seq->seq[length++] = (char)c;
	    }

	}
	if (c == '>') ungetc(c,fp);
	
	seq->seq[length] = '\0';
	//seq->max         = max;
	seq->length      = length;
	
	if (length > 3000) //ZAM TODO remove magic number 3000 
	  {
	    printf("We have met a read which is longer than our max expected, of 3000 bp. Legnth is %d, Exiting.", length);
	    exit(1);
	  }

	return 1;
      }
    else{
      printf("Not fasta");
      exit(1);
      return 0;
    }
  }

  return 0;
}



//caller of this takes ownership of the Seqeucne* that is malloc-ed in ere
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

    if (line[0] == '@')
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
	
	
	while (!feof(fp) && ((c = fgetc(fp)) != '\n')   ) {
	  if (isalpha(c)) {
	    if (length + 1 >= max) {
	      max += SEQ_BLOCK_SIZE;
	      seq->seq = (char*)realloc(seq->seq, sizeof(char) * max);
	      if (seq->seq == NULL)
		{
		  printf("Out of mem doing realloc of seq");
		  exit(1);
		}

	    }
	  seq->seq[length++] = (char)c;
	  //printf("Next seq character is %c and length is now %d\n",c, length);
	  }
	  else
	    {
	      printf("invalid fastq");
	      exit(1);
	    }
	}
	
	//ignore strand
	char temp[3];
	fgets(temp,3,fp);
	//printf("strand %sZam",temp);
	
	//get quality - we know need same amount as seq->seq
	int qual_length=0;
	seq->qual = (char*)realloc(seq->qual, sizeof(char) * max);
	if (seq->qual == NULL)
	  {
	    printf("Out of mem doing realloc of qual");
	    exit(1);
	  }

	
	while ( !feof(fp) && ((c = fgetc(fp)) != '\n') ) {
	seq->qual[qual_length++] = (char)c;
	//printf("Next qual character is %c and qual_length is now %d\n",c, qual_length);
	}
	
	//printf("After getting quality, qual length is %d\n", qual_length);
	seq->seq[length] = '\0';
	//printf("Sequence is %s,length is %d\n",seq->seq, length);
	seq->qual[qual_length]='\0';
	//printf("Quality %sZam",seq->qual);
	seq->max         = max;
	seq->length      = length;
	
	
	
	if ((qual_length != length))
	  {
	    printf("Not reading fastq correctly, qual length is %d and seq length is %d\nqual is %s and sequence is %s", qual_length, length, seq->qual, seq->seq);
	    exit(1);
	  }
	
	
	
	return seq;
      }
    else{
      free(seq);
      return NULL;
    }
  }
  
  free(seq);
  return NULL;
}




void free_sequence(Sequence ** sequence){
  free((*sequence)->seq);
  free((*sequence)->qual);
  free(*sequence);
  *sequence = NULL;
}

