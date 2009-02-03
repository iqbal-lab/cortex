#include <stdio.h>
#include <stdlib.h>
#include <seq.h>
#include <limits.h>


/* Read sequence from file "fp" in FASTA format.
   it returns the length of the sequence, 0 if no sequence is left in file
*/

int read_sequence_from_fasta(FILE *fp, Sequence * seq, int max_read_length){

  char line[LINE_MAX]; //LINE_MAX is defined in limits.h
  int i;
  int j = 0; //length of sequence
  long file_pointer;

  if (fp == NULL){
    fputs("File not defined\n",stderr);
    exit(1);
  }

  if (seq == NULL){
    fputs("Cannot pass a NULL pointer for seq\n",stderr);
    exit(1);
  }

  if (seq->seq == NULL){
    fputs("Dont give me a null pointer for seq->seq - alloc memory yourself and give me that\n",stderr);
    exit(1);
  }

  if (seq->qual == NULL){
    fputs("Dont give me a null pointer for seq->qual - alloc memory yourself and give me that\n",stderr);
    exit(1);
  }

  if (seq->name == NULL){
    fputs("Dont give me a null pointer for seq->name - alloc memory yourself and give me that\n",stderr);
    exit(1);
  }

  
  //read name of fasta entry
  if (fgets(line, LINE_MAX, fp) != NULL){
    if (line[0] == '>'){
      for(i = 1;i<LINE_MAX;i++){	   
	if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r'){
	  break;
	}
	if(i>255){
	  fputs("Name too long\n",stderr);
	  exit(1);
	}
	seq->name[i-1] = line[i];
      }
      seq->name[i-1] = '\0';
    
      //read sequence -- verify position first
      file_pointer = ftell(fp);

      while (fgets(line, LINE_MAX, fp) != NULL){
	if (line[0] == '>'){
	  fseek(fp,file_pointer,SEEK_SET);
	  break;
	}
	
	//check line is fine
	for(i=0;i<LINE_MAX;i++){
	  if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r'){
	    break; //fine but nothing to add
	  }
	  
	  if (line[i] != 'A' && line[i] != 'a' && 
	      line[i] != 'C' && line[i] != 'c' && 
	      line[i] != 'G' && line[i] != 'g' && 
	      line[i] != 'T' && line[i] != 't' && 
	      line[i] != 'N' && line[i] != 'n'){ 
	    
	    fprintf(stderr,"Invalid symbol [%c] in entry\n",line[i]);
	    exit(1);
	  }
	  
	  seq->seq[j]  = line[i];
	  seq->qual[j] = '\0';
	  j++;
	  
	  if (j==max_read_length){
	    fprintf(stderr,"read [%s] too long [%i]\n",seq->name,j);
	    exit(1);
	  }
	  
	}
    
	file_pointer = ftell(fp);
      }
    }
    else{
      fputs("syntax error in fasta file -- it misses >\n",stderr);
      exit(1);
    }
  }
  seq->seq[j]  = '\0';
  seq->qual[j] = '\0';
  return j;
}



/* Read sequence from file "fp" in FASTQ format.
   it returns the length of the sequence, 0 if no sequence is left in file
*/

int read_sequence_from_fastq(FILE *fp, Sequence * seq, int max_read_length){

  char line[LINE_MAX]; //LINE_MAX is defined in limits.h
  int i;
  int j = 0; //length of sequence
  int q = 0; //length of qualities
  long file_pointer;

  if (fp == NULL){
    fputs("File not defined\n",stderr);
    exit(1);
  }

  if (seq == NULL){
    fputs("Cannot pass a NULL pointer for seq\n",stderr);
    exit(1);
  }

  if (seq->seq == NULL){
    fputs("Dont give me a null pointer for seq->seq - alloc memory yourself and give me that\n",stderr);
    exit(1);
  }

  if (seq->name == NULL){
    fputs("Dont give me a null pointer for seq->name - alloc memory yourself and give me that\n",stderr);
    exit(1);
  }

  if (seq->qual == NULL){
    fputs("Dont give me a null pointer for seq->qual - alloc memory yourself and give me that\n",stderr);
    exit(1);
  }

  
  //read name of fastq entry
  if (fgets(line, LINE_MAX, fp) != NULL){
    if (line[0] == '@'){
      for(i = 1;i<LINE_MAX;i++){	   
	if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r'){
	  break;
	}
	if(i>255){
	  fputs("Name too long\n",stderr);
	  exit(1);
	}
	seq->name[i-1] = line[i];
      }
      seq->name[i-1] = '\0';
    
      //read sequence 

      while (fgets(line, LINE_MAX, fp) != NULL){

	if ((line[0] == '+') || (line[0] == '-')){ //go to get qualities
	  break;
	}
	
	//check line is fine
	for(i=0;i<LINE_MAX;i++){
	  if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r'){
	    break; //fine but nothing to add
	  }
	  
	  if (line[i] != 'A' && line[i] != 'a' && 
	      line[i] != 'C' && line[i] != 'c' && 
	      line[i] != 'G' && line[i] != 'g' && 
	      line[i] != 'T' && line[i] != 't' && 
	      line[i] != 'N' && line[i] != 'n'){ 
	    
	    fprintf(stderr,"Invalid symbol [%c] in entry\n",line[i]);
	    exit(1);
	  }
	  
	  seq->seq[j] = line[i];
	  j++;
	  
	  if (j==max_read_length){
	    fprintf(stderr,"read [%s] too long [%i]\n",seq->name,j);
	    exit(1);
	  }
	  
	}
      }

      //read qualities -- verify position first
      file_pointer = ftell(fp);

      while (fgets(line, LINE_MAX, fp) != NULL){
	if (line[0] == '@' && j==q){
	  fseek(fp,file_pointer,SEEK_SET);
	  break;
	}
       
	for(i=0;i<LINE_MAX;i++){
	  if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r'){
	    break; //fine but nothing to add
	  }
	  
	  seq->qual[q] = line[i];
	  q++;
	  
	  if (q==max_read_length){
	    fprintf(stderr,"qualities for [%s] too long [%i]\n",seq->name,q);
	    exit(1);
	  }
	  
	}
    
	file_pointer = ftell(fp);
      }

      if (j!=q){
	fprintf(stderr,"qualities [%i] and sequence [%i] sizes don't coincide for [%s]\n",q,j,seq->name);
	exit(1);
      }

    }
    else{
      fputs("syntax error in fastq file -- it misses @\n",stderr);
      exit(1);
    }
  }
  seq->seq[j]   = '\0';
  seq->qual[q]  = '\0'; //this is not technical necessary but simplyfies many checks downstream
  return j;
}

void free_sequence(Sequence ** sequence){
  free((*sequence)->seq);
  if ((*sequence)->qual != NULL){
    free((*sequence)->qual);
  }
  free(*sequence);
  *sequence = NULL;
}

