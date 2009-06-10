#include <stdio.h>
#include <stdlib.h>
#include <seq.h>
#include <limits.h>
#include <string.h>

/* Read sequence from file "fp" in FASTA format.
   it returns the length of the sequence, 0 if no sequence is left in file
*/
boolean good_base(char c);


char current_entry_name[LINE_MAX+1] = "";
int  last_end_coord;

//backtrack updates the file pointer to adjust for next kmer, eg
//>pepe
//AACCCTGCTGC
//first time you call 
// read_sequence_from_long_fasta(fp, seq, 3, &full_entry,1)
//    reads AAC and bractracks file pointer to start at C 
// then read_sequence_from_long_fasta(fp, seq, 4, &full_entry,0)
//returns CCCT

//new_entry tells the parser to expect a new fasta entry (ie starts with >)
//full_entry tells the caller if the end of the entry has been reached

int read_sequence_from_fasta(FILE *fp, Sequence * seq, int max_chunk_length,boolean new_entry, boolean * full_entry, int offset){

  char line[LINE_MAX]; //LINE_MAX is defined in limits.h
  int i;
  int j = 0; //length of sequence
  boolean good_read;

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

 
  //get name
  if (new_entry == true){ //need a '>' follwed by a name
    if (fgets(line, LINE_MAX, fp) != NULL){
      if (line[0] == '>'){
	i=1;
	while(i<LINE_MAX && ! (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r')){
	  seq->name[i-1] = line[i];
	  i++;
	}
	seq->name[i-1] = '\0';
      }
      else{
	fprintf(stderr,"syntax error in fasta entry %s\n",line);
	exit(1);
      }
    }
    

    strcpy(current_entry_name,seq->name);
    seq->start = 1; //new entry
  }
  else{ //old entry
    strcpy(seq->name,current_entry_name);
    seq->start = last_end_coord+1;
  }

  //printf("entry: %s %i\n",seq->name,seq->start);
    
  //get sequence

  boolean ready_to_return = false;
  long file_pointer = 0;
  long prev_file_pointer = 0;

  j=offset; //global length
  *full_entry=false;

  file_pointer = ftell(fp);
  

  while (ready_to_return==false){

    file_pointer = ftell(fp);

    if (fgets(line,LINE_MAX,fp) != NULL){

      int length = strlen(line);

      //sanity check
      
      if (line[0] == '>'){
	fseek(fp,file_pointer,SEEK_SET);
	ready_to_return = true;
	*full_entry = true;
	//printf("complete line - new entry\n");
      }
      else{
	if (j==max_chunk_length){
	  fseek(fp,prev_file_pointer,SEEK_SET); 
	  ready_to_return = true;
	  *full_entry = false;
	  //printf("complete line - same entry\n");
	}
	else{
	  //check line is fine
	  i=0; //counter within line
	  while(i<length && ! (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r') && ready_to_return==false){
	    good_read = good_base(line[i]);	  
	    if (! good_read){
	      fprintf(stderr,"Invalid symbol [%c] pos:%i in entry %s\n",line[i],i,seq->name);
	      //exit(1);
	    }
	    
	    seq->seq[j]  = line[i];
	    seq->qual[j] = '\0';
	    i++;
	    file_pointer++;
	    j++;
	    
	    if (j==max_chunk_length){	    
	      //check if line is not complete
	      if (i<length-1){
		ready_to_return = true;
		*full_entry = false;
		fseek(fp,file_pointer,SEEK_SET); 	      
		//printf("incomplete line - incomplete entry\n");
	      }

	      //else 3 cases might happen with line complete
	      // 1. next line starts ">" -> *full_entry = true
	      // 2. next line is EOF -> *full_entry = true
	      // 3. next line is more sequence -> *full_entry = false
	    }
	  }	    
	}	    
      }
    }
    else{
      *full_entry = true;
      ready_to_return = true;
      //printf("complete line - complete entry - complete file\n");
    }
    
    prev_file_pointer = file_pointer;
  }
    
  seq->end = seq->start+j-1-offset;
  last_end_coord = seq->end;
  
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
  boolean good_read = true;
  

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

  do{
    //read name of fastq entry
    if (fgets(line, LINE_MAX, fp) != NULL){
      if (line[0] == '@'){
	for(i = 1;i<LINE_MAX;i++){	   
	  if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r'){
	    break;
	  }
	  if(i>LINE_MAX){
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
	    
	    good_read = good_base(line[i]);
	    if(! good_read){
	      fprintf(stderr,"Invalid symbol [%c]  pos:%i in entry %s\n",line[i],i,seq->name);
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
  } while (! good_read); 

  seq->seq[j]   = '\0';
  seq->qual[q]  = '\0'; //this is not technical necessary but simplyfies many checks downstream
  return j;
}


void alloc_sequence(Sequence * seq, int max_read_length, int max_name_length){
 
 
  if (seq == NULL){							
    fputs("Cannot pass a null seq to sequence_alloc\n",stderr);	
    exit(1);								
  }
  seq->name = malloc(sizeof(char) * max_name_length);
  if (seq->name == NULL){
    fputs("Out of memory trying to allocate string\n",stderr);
    exit(1);
  }
  seq->seq  = malloc(sizeof(char) * (max_read_length+1));
  if (seq->seq == NULL){
    fputs("Out of memory trying to allocate string\n",stderr);
    exit(1);
  }
  seq->qual  = malloc(sizeof(char) * (max_read_length+1));
  if (seq->qual == NULL){
    fputs("Out of memory trying to allocate string\n",stderr);
    exit(1);
  }
}


boolean good_base(char c){
  boolean ret;
  if (c  != 'A' && c != 'a' && 
      c != 'C' && c != 'c' && 
      c != 'G' && c != 'g' && 
      c != 'T' && c != 't' && 
      c != 'N' && c != 'n' 
      ){
    ret = false;
  }	
  else{
    ret =  true;
  }
  
  return ret;
}


void free_sequence(Sequence ** sequence){
  free((*sequence)->name);
  free((*sequence)->seq);
  free((*sequence)->qual);
  free(*sequence);
  *sequence = NULL;
}

void shift_last_kmer_to_start_of_sequence(Sequence * sequence, int length, short kmer_size){

  int i;

  if (length-kmer_size<kmer_size){
    puts("kmer_size too long\n");
    exit(1);
  }

  for(i=0;i<kmer_size; i++){
    sequence->seq[i] = sequence->seq[length-kmer_size+i];   
    sequence->qual[i] = '\0';
  }

}
