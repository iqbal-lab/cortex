
#include <stdio.h>
#include <stdlib.h>
#include <seq.h>
#include <limits.h>


/* Read sequence from file "fp" in FASTA format.
   it returns the length of the sequence, 0 if no sequence is left in file
*/
boolean good_base(char c);

int read_sequence_from_fasta(FILE *fp, Sequence * seq, int max_read_length){

  char line[LINE_MAX]; //LINE_MAX is defined in limits.h
  int i;
  int j = 0; //length of sequence
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

  if (seq->qual == NULL){
    fputs("Dont give me a null pointer for seq->qual - alloc memory yourself and give me that\n",stderr);
    exit(1);
  }

  if (seq->name == NULL){
    fputs("Dont give me a null pointer for seq->name - alloc memory yourself and give me that\n",stderr);
    exit(1);
  }

  boolean end_of_file=false;

  do{
    //read name of fasta entry
    if (fgets(line, LINE_MAX, fp) != NULL)
      {
	if (line[0] == '>')
	  {
	    for(i = 1;i<LINE_MAX;i++)
	      {	   
		if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r')
		  {
		    break;
		  }
		if(i>LINE_MAX)
		  {
		    fputs("Name too long\n",stderr);
		    exit(1);
		  }
		seq->name[i-1] = line[i];
	      }
	    seq->name[i-1] = '\0';
	    
	    //read sequence -- verify position first
	    file_pointer = ftell(fp);
	    
	    while (fgets(line, LINE_MAX, fp) != NULL)
	      {
		if (line[0] == '>')
		  {
		    fseek(fp,file_pointer,SEEK_SET);
		    break;
		  }
		
		//check line is fine
		for(i=0;i<LINE_MAX;i++)
		  {
		    if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r')
		      {
			break; //fine but nothing to add
		      }
		    
		    good_read = good_base(line[i]);
		    
		    if (! good_read)
		      {
			//fprintf(stdout,"Invalid symbol [%c] pos:%i in entry %s. Ignoring and moving on to next read.\n",line[i],i,seq->name);
			good_read=false;
			j=0;
			break;		  //exit(1);
		      }
		    
		    seq->seq[j]  = line[i];
		    seq->qual[j] = '\0';
		    j++;
		    
		    if (j==max_read_length)
		      {
			//fprintf(stdout,"read [%s] too long [%i]. Ignore and move on to next\n",seq->name,j);
			// exit(1);
			good_read=false;
			j=0;
			break;
		      }
		    
		  }
		
		file_pointer = ftell(fp);
	      }
	    
	    
	  }
	else
	  {
	    fputs("syntax error in fasta file -- it misses >\n",stderr);
	    exit(1);
	  }
      }
    else
      {
	end_of_file=true;
      }
    
  } while (! good_read && !end_of_file);
  
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

  boolean end_of_file=false;

  do{

    good_read=true;
    q=0;
    j=0;

    //read name of fastq entry
    if (fgets(line, LINE_MAX, fp) != NULL)
      {
	if (line[0] == '@')
	  {
	    for(i = 1;i<LINE_MAX;i++)
	      {	   
		if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r')
		  {
		    break;
		  }
		if(i>LINE_MAX)
		  {
		    fputs("Name too long\n",stderr);
		    exit(1);
		  }
		seq->name[i-1] = line[i];
	      }
	    seq->name[i-1] = '\0';
	    
	    //read sequence 
	    
	    while (fgets(line, LINE_MAX, fp) != NULL)
	      {
		
		if ((line[0] == '+') || (line[0] == '-'))
		  { //go to get qualities
		    break;
		  }
		
		//check line is fine
		for(i=0;i<LINE_MAX;i++)  
		  {
		    if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r')
		      {
			break; //fine but nothing to add
		      }
		    
		    if (! good_base(line[i]))
		    {
		      good_read=false;
		      //fprintf(stdout,"Invalid symbol [%c]  pos:%i in entry %s. Will ignore this read.\n",line[i],i,seq->name);
		    }


		    seq->seq[j] = line[i];
		    j++;
		    
		    if (j==max_read_length){
		      fprintf(stdout,"read [%s] too long [%i]. Exiting...\n",seq->name,j);
		      exit(1);
		    }
		    
		  }
	      }
	    
	    //read qualities -- verify position first
	    file_pointer = ftell(fp);
	    
		
	    while (fgets(line, LINE_MAX, fp) != NULL)
	      {

		if (line[0] == '@' && (j<=q) )//then we have gone on to the next read  
		  //allowing q>j in case where qualities longer than j
		  {
		    fseek(fp,file_pointer,SEEK_SET);
		    break; //goto next read
		  }
		
		
		for(i=0;i<LINE_MAX;i++)
		  {
		    if (line[i] == '\n' || line[i] == ' ' || line[i] == '\t' || line[i] == '\r'){
		      break; //fine but nothing to add
		    }
		    
		    seq->qual[q] = line[i];
		    q++;
		    
		    if (q==max_read_length)
		      {
			fprintf(stdout,"qualities for [%s] longer than the max read length  [%i]. Exiting...\n",seq->name,q);
			exit(1);
		      }
		    
		  }
		
		file_pointer = ftell(fp);
	      }
	    
	    
	    if (j!=q)
	      {
		fprintf(stdout,"qualities [%i] and sequence [%i] sizes don't coincide for [%s]. Skip it\n",q,j,seq->name);
		good_read=false;
	      }
	    
	  }//if line starts with @
	else
	  {
	    fputs("syntax error in fastq file -- it misses @\n",stderr);
	    exit(1);
	  }
      }
    else
      {
	end_of_file=true;
      }

  } while (! good_read && !end_of_file); 

  seq->seq[j]   = '\0';
  seq->qual[q]  = '\0'; //this is not technically necessary but simplifies many checks downstream
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
      c != 'N' && c != 'n'){       
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

