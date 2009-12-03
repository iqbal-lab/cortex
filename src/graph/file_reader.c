/*
  routines to load files into dB Graph
 */



#include <binary_kmer.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <dB_graph.h>
#include <seq.h>
#include <file_reader.h>
#include <element.h>
#include <string.h>

int MAX_READ_LENGTH=1000;
int MAX_FILENAME_LENGTH=200;



long long load_seq_into_graph(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
			      long long * bad_reads, char qualiy_cut_off, long long* dup_reads,
			      int max_read_length, boolean remove_dups_single_endedly, dBGraph * db_graph);

long long load_paired_end_seq_into_graph(FILE* fp1, FILE* fp2, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry),
                                         long long * bad_reads, char quality_cut_off, int max_read_length, long long* dup_reads, boolean remove_dups, dBGraph * db_graph);



long long load_fastq_from_filename_into_graph(char* filename, long long * bad_reads,  char quality_cut_off, long long* dup_reads, int max_read_length, 
					      boolean remove_duplicates_single_endedly, dBGraph* db_graph)
{

  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      puts("new_entry has to be true for fastq\n");
      exit(1);
    }

    return read_sequence_from_fastq(fp,seq,max_read_length);
  }

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    fprintf(stderr,"cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  long long ret;

  if (remove_duplicates_single_endedly==true)
    {
      ret =  load_seq_into_graph(fp,&file_reader,bad_reads,quality_cut_off,dup_reads, max_read_length,true, db_graph);
    }
  else
    {
      ret =  load_seq_into_graph(fp,&file_reader,bad_reads,quality_cut_off, dup_reads, max_read_length,false,db_graph);
    }
  fclose(fp);

  return ret;
}

long long load_paired_fastq_from_filenames_into_graph(char* filename1, char* filename2, 
						      long long * bad_reads,  char quality_cut_off, int max_read_length, 
						      long long* dup_reads, boolean remove_duplicates, dBGraph* db_graph)
{

  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      puts("new_entry has to be true for fastq\n");
      exit(1);
    }

    return read_sequence_from_fastq(fp,seq,max_read_length);
  }

  FILE* fp1 = fopen(filename1, "r");
  if (fp1 == NULL){
    fprintf(stderr,"cannot open file:%s\n",filename1);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }
  FILE* fp2 = fopen(filename2, "r");
  if (fp2 == NULL){
    fprintf(stderr,"cannot open file:%s\n",filename2);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  long long ret =  load_paired_end_seq_into_graph(fp1, fp2, &file_reader,bad_reads,quality_cut_off,max_read_length, dup_reads, remove_duplicates, db_graph);
  fclose(fp1);
  fclose(fp2);

  //printf("Total duplicate reads in %s and %s: %lld\n", filename1, filename2, *dup_reads);
  return ret;
}


//this routine supports big fasta entries (chromosome length for example)
long long load_fasta_from_filename_into_graph(char* filename, long long * bad_reads, long long* dup_reads, int max_chunk_length, 
					      boolean remove_duplicates_single_endedly, dBGraph* db_graph)
{

  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);

    return ret;
  }

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    fprintf(stderr,"cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  long long ret;
  if (remove_duplicates_single_endedly==true)
    {
      ret =  load_seq_into_graph(fp,&file_reader,bad_reads,0,dup_reads, max_chunk_length,true, db_graph);
    }
  else
    {
      ret =  load_seq_into_graph(fp,&file_reader,bad_reads,0,dup_reads, max_chunk_length,false, db_graph);
    }
  fclose(fp);


  return ret;
}



//do not export the folloiwing internal function
void load_kmers_from_sliding_window_into_graph_marking_read_starts(KmerSlidingWindowSet * windows, boolean* prev_full_ent, 
								   boolean* full_ent, long long* seq_len, boolean mark_read_starts, dBGraph* db_graph)
{

      Element * current_node  = NULL;
      Element * previous_node  = NULL;
      Orientation current_orientation=forward;
      Orientation previous_orientation=forward;
      BinaryKmer tmp_kmer;
      int i,j;
      for(i=0;i<windows->nwindows;i++){ //for each window
	KmerSlidingWindow * current_window = &(windows->window[i]);
	
	for(j=0;j<current_window->nkmers;j++){ //for each kmer in window
	  boolean found = false;
	  current_node = hash_table_find_or_insert(element_get_key(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),&found,db_graph);	  
	  if (current_node == NULL){
	    fputs("file_reader: problem - current kmer not found\n",stderr);
	    exit(1);
	  }
	  
	  if (! (i==0 && j==0 && *prev_full_ent == false && current_node == previous_node)){ //otherwise is the same old last entry
	    element_update_coverage(current_node,1);
	  }

	  current_orientation = db_node_get_orientation(&(current_window->kmer[j]),current_node, db_graph->kmer_size);


	  if (mark_read_starts==true)
	    {
	      if ( (i==0) && (j==0) && (current_orientation==forward))
		{
		  db_node_set_read_start_status(current_node, forward);
		}
	      else if ( (i==0) && (j==0) && (current_orientation==reverse))
		{
		  db_node_set_read_start_status(current_node, reverse);
		}
	    }

	  
	  if (DEBUG){
	    char kmer_seq[db_graph->kmer_size];
	    char kmer_seq2[db_graph->kmer_size];
	    
	    printf("kmer i:%i j:%i:  %s %s %i\n",i,j,binary_kmer_to_seq(&(current_window->kmer[j]),db_graph->kmer_size,kmer_seq),
		   binary_kmer_to_seq(binary_kmer_reverse_complement(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph->kmer_size,kmer_seq2),element_get_coverage(current_node));
	  }
	  
	  if (j>0){
	    
	    if (previous_node == NULL){
	      printf("PROBLEM i:%i j:%i seq_length:%qd nkmers:%i prev_full_entry:%s\n",i,j,*(seq_len),current_window->nkmers,*prev_full_ent == true ? "true" : "false");
	      fputs("file_reader: problem - prev kmer not found\n",stderr);
	      exit(1);
	    }
	    else{
	      db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size);	  	  
	    }
	  }
	  previous_node = current_node;
	  previous_orientation = current_orientation;
	  
	}
      }
     
}



//returns length of sequence loaded
//penultimate argument is to specify whether to discard potential PCR duplicate reads single-endedly - ie if read starts at same kmer
// as a previous read, then discard it. This is a pretty harsh filter, and if ppssible, prefer to use paired end info.
// So in general when calling this function, would expect that boolean to be set to false, unless you know you have low coverage, so have
// low probability of two reads starting at the same point.
long long load_seq_into_graph(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
			      long long * bad_reads, char quality_cut_off, long long* dup_reads, 
			      int max_read_length, boolean remove_dups_single_endedly, dBGraph * db_graph){

 
  //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  
  long long seq_length=0;
  short kmer_size = db_graph->kmer_size;

  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);
  
  int entry_length;
  boolean full_entry = true;
  boolean prev_full_entry = true;
  Element * previous_node = NULL;

  
  while ((entry_length = file_reader(fp,seq,max_read_length,full_entry,&full_entry))){
   
    //printf("entry %s %i %i %s\n",seq->name,seq->start,seq->end,prev_full_entry?"not connect":"connect");
    
 
   if (DEBUG){
      printf ("\nsequence %s - kmer size: %i - entry length: %i - max kmers:%i\n",seq->seq,db_graph->kmer_size, entry_length, max_kmers);
    }
    
    int i,j;
    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));
   
    //printf("Length %qd %i %s %i %i\n",seq_length,entry_length, seq->name, seq->start, seq->end);

 
    int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,
						   entry_length,quality_cut_off,db_graph->kmer_size,windows,max_windows, max_kmers);

    //printf("number kmers:%i %i\n",nkmers,windows->nwindows);


    if (remove_dups_single_endedly==true)
      {
	BinaryKmer tmp_kmer;
	KmerSlidingWindow * first_window= &(windows->window[0]);
	//get node and orientation for first kmer in read from file1
	dBNode* test = hash_table_find(element_get_key(&(first_window->kmer[0]),db_graph->kmer_size, &tmp_kmer), db_graph);

	if (test != NULL) // can only be previous read_start if the node already s there in the graph!
	  {
	    Orientation test_o = db_node_get_orientation(&(first_window->kmer[0]),test, db_graph->kmer_size);

	    if (db_node_check_single_ended_duplicates(test, test_o)==true)
	      {
		(*dup_reads)++;
		//printf("Discard reads %s and %s - duplicates\n", seq1->name, seq2->name);
		continue;
	      }
	  }
      }

    
    if (nkmers == 0) {
      (*bad_reads)++;
    }
    else {
      if (remove_dups_single_endedly==true)
	{
	  load_kmers_from_sliding_window_into_graph_marking_read_starts(windows, &prev_full_entry, &full_entry, &seq_length, true, db_graph);
	}
      else
	{
	  //don't mark read-starts
	  load_kmers_from_sliding_window_into_graph_marking_read_starts(windows, &prev_full_entry, &full_entry, &seq_length, false, db_graph);
	}
    }

    if (full_entry == false){
      shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
    }

    prev_full_entry = full_entry;
   

  }
  
 
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
  return seq_length;    
}



//do not export the following internal function
void paired_end_sequence_core_loading_loop(FILE* fp1 , FILE* fp2, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry),
					   Sequence* seq1, Sequence* seq2, char quality_cut_off, int max_read_length, int max_kmers, int max_windows, 
					   KmerSlidingWindowSet * windows1, KmerSlidingWindowSet * windows2, long long* seq_len, long long* dup_reads, long long* bad_reads, 
					   boolean remove_dups, dBGraph* db_graph)
{

  int entry_length1=0;
  int entry_length2=0;

  boolean full_entry1 = true;
  boolean full_entry2 = true;
  boolean prev_full_entry1 = true;
  boolean prev_full_entry2 = true;

 
  //we assume that both files are the same length
  while ((entry_length1 = file_reader(fp1,seq1,max_read_length,full_entry1,&full_entry1))){
   
    entry_length2 = file_reader(fp2,seq2,max_read_length,full_entry2,&full_entry2);

    if (entry_length2==0)
      {
	printf("Warning - second mate pair file is shorter than first. Last read was %s\n", seq1->name);
	exit(1);
      }
    

    //printf("entry %s %i %i %s\n",seq1->name,seq1->start,seq1->end,prev_full_entry1?"not connect":"connect");
    //printf("entry %s %i %i %s\n",seq2->name,seq2->start,seq2->end,prev_full_entry2?"not connect":"connect");
    
   if (DEBUG){
      printf ("\nsequence 1 of pair %s - kmer size: %i - entry length: %i - max kmers:%i\n",seq1->seq,db_graph->kmer_size, entry_length1, max_kmers);
      printf ("\nsequence 2 of pair %s - kmer size: %i - entry length: %i - max kmers:%i\n",seq2->seq,db_graph->kmer_size, entry_length2, max_kmers);
    }
    
    //seq_len is ampount fo sequence read in from files, not necessarily loaded into graph
    (*seq_len) += (long long) (entry_length1 - (prev_full_entry1==false ? db_graph->kmer_size : 0));
    (*seq_len) += (long long) (entry_length2 - (prev_full_entry2==false ? db_graph->kmer_size : 0));
   
 
    int nkmers1 = get_sliding_windows_from_sequence(seq1->seq,seq1->qual,
						   entry_length1,quality_cut_off,db_graph->kmer_size,windows1,max_windows, max_kmers);
    int nkmers2 = get_sliding_windows_from_sequence(seq2->seq,seq2->qual,
						   entry_length2,quality_cut_off,db_graph->kmer_size,windows2,max_windows, max_kmers);


    //check whether first kmer of both reads is a read-start. If yes, then discard read as duplicate.
    BinaryKmer tmp_kmer;

    if (remove_dups==true)
      {
	KmerSlidingWindow * first_window1= &(windows1->window[0]);
	KmerSlidingWindow * first_window2= &(windows2->window[0]);
	//get node and orientation for first kmer in read from file1
	dBNode* test1 = hash_table_find(element_get_key(&(first_window1->kmer[0]),db_graph->kmer_size, &tmp_kmer), db_graph);
	//get node and orientation for first kmer in read from file2
	dBNode* test2 = hash_table_find(element_get_key(&(first_window2->kmer[0]),db_graph->kmer_size, &tmp_kmer), db_graph);
	

	if ((test1 != NULL) && (test2 !=NULL))// can only be previous read_starts if the node already s there in the graph!
	  {
	    Orientation test_o1 = db_node_get_orientation(&(first_window1->kmer[0]),test1, db_graph->kmer_size);
	    Orientation test_o2 = db_node_get_orientation(&(first_window2->kmer[0]),test2, db_graph->kmer_size);

	    if (db_node_check_duplicates(test1, test_o1, test2, test_o2)==true)
	      {
		(*dup_reads)++;
		//printf("Discard reads %s and %s - duplicates\n", seq1->name, seq2->name);
		continue;
	      }
	  }
      }
    

    if ((nkmers1 == 0)&& (nkmers2 == 0)) {
      (*bad_reads)=(*bad_reads)+2;
    }
    else if ((nkmers1==0)&&(nkmers2!=0))
      {
	(*bad_reads)++;
	load_kmers_from_sliding_window_into_graph_marking_read_starts(windows2, &prev_full_entry2, &full_entry2, seq_len, true, db_graph);
      }
    else if ((nkmers1!=0)&&(nkmers2==0))
      {
	(*bad_reads)++;
	load_kmers_from_sliding_window_into_graph_marking_read_starts(windows1, &prev_full_entry1, &full_entry1, seq_len, true, db_graph);
      }
    else 
      {
	load_kmers_from_sliding_window_into_graph_marking_read_starts(windows1, &prev_full_entry1, &full_entry1, seq_len, true, db_graph);
	load_kmers_from_sliding_window_into_graph_marking_read_starts(windows2, &prev_full_entry2, &full_entry2, seq_len, true, db_graph);
      }

    }

    if (full_entry1 == false){
      shift_last_kmer_to_start_of_sequence(seq1,entry_length1,db_graph->kmer_size);
    }
    if (full_entry2 == false){
      shift_last_kmer_to_start_of_sequence(seq2,entry_length2,db_graph->kmer_size);
    }


    prev_full_entry1 = full_entry1;
    prev_full_entry2 = full_entry2;
   

}
  


//returns length of sequence loaded
long long load_paired_end_seq_into_graph(FILE* fp1, FILE* fp2, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
					 long long * bad_reads, char quality_cut_off, int max_read_length, long long* dup_reads, boolean remove_dups, dBGraph * db_graph){

 


  //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq1 = malloc(sizeof(Sequence));
  if (seq1 == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq1,max_read_length,LINE_MAX);
  Sequence * seq2 = malloc(sizeof(Sequence));
  if (seq2 == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq2,max_read_length,LINE_MAX);
  
  
  long long seq_length=0;
  short kmer_size = db_graph->kmer_size;

  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows1 = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows1 == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows1, max_windows, max_kmers);

  KmerSlidingWindowSet * windows2 = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows2 == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows2, max_windows, max_kmers);
  

  //work through the paired files, loading sequence into the graph. If remove_dups==true, then for each pair, discard both mate reads if they both start at a read_start kmer
  paired_end_sequence_core_loading_loop(fp1, fp2, file_reader, seq1, seq2, quality_cut_off,  max_read_length, max_kmers, max_windows, 
					windows1, windows2,  &(seq_length), dup_reads, bad_reads, remove_dups, db_graph);
  


  free_sequence(&seq1);
  free_sequence(&seq2);
  binary_kmer_free_kmers_set(&windows1);
  binary_kmer_free_kmers_set(&windows2);
  return seq_length;    
}



//assume we have two lists of equal length, of mate fastq files in the same order in each file
long long load_list_of_paired_end_fastq_into_graph(char* list_of_left_mates, char* list_of_right_mates, char quality_cut_off, int max_read_length, 
						   long long* bad_reads, long long* num_dups, boolean remove_dups, dBGraph* db_graph)
{


  FILE* f1 = fopen(list_of_left_mates, "r");
  if (f1==NULL)
    {
      printf("Cannot open %s - exit\n", list_of_left_mates);
      exit(1);
    }
  FILE* f2 = fopen(list_of_right_mates, "r");
  if (f2==NULL)
    {
      printf("Cannot open %s - exit\n", list_of_right_mates);
      exit(1);
    }

  long long seq_length = 0;
  char filename1[MAX_FILENAME_LENGTH+1];
  char filename2[MAX_FILENAME_LENGTH+1];


  while (fgets(filename1,MAX_FILENAME_LENGTH, f1) !=NULL);
    {
      //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(filename1, '\n')) != NULL)
	*p = '\0';

      if (fgets(filename2,MAX_FILENAME_LENGTH, f2) ==NULL)
	{
	  printf("Files %s and %s are not the same length, and they are lists of mate fastqs\n",list_of_left_mates, list_of_right_mates);
	}
      if ((p = strchr(filename2, '\n')) != NULL)
	*p = '\0';


      seq_length += load_paired_fastq_from_filenames_into_graph(filename1, filename2, bad_reads, quality_cut_off, max_read_length, num_dups, remove_dups, db_graph);

    }



    return seq_length;

}

//returns number of sequence loaded (ie all the kmers concatenated)
//count_kmers returns the number of new kmers

long long load_binary_from_filename_into_graph(char* filename,  dBGraph* db_graph){

  FILE* fp_bin = fopen(filename, "r");
 
  dBNode node_from_file;
  boolean found;
  long long count=0;
  BinaryKmer tmp_kmer;

  if (fp_bin == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }
  
  //Go through all the entries in the binary file
  while (db_node_read_binary(fp_bin,db_graph->kmer_size,&node_from_file)){
    count++;
    
    dBNode * current_node  = NULL;
    current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size, &tmp_kmer),&found,db_graph);
       
    db_node_set_edges(current_node,db_node_get_edges(&node_from_file));
    element_update_coverage(current_node, element_get_coverage(&node_from_file));
  }
 
  fclose(fp_bin);

  return (long long) db_graph->kmer_size * count;

}


void read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
									    long long * bad_reads, char quality_cut_off, int max_read_length, dBGraph * db_graph)
{
    //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  
  int seq_length=0;
  short kmer_size = db_graph->kmer_size;

  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);

  BinaryKmer tmp_kmer;

  boolean full_entry = true;
  boolean prev_full_entry = true;
  
  int entry_length;

  while ((entry_length = file_reader(fp,seq,max_read_length, full_entry, &full_entry))){

    if (DEBUG){
      printf ("\nsequence %s\n",seq->seq);
    }
    
    int i,j;
    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));

    
    int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,quality_cut_off,db_graph->kmer_size,windows,max_windows, max_kmers);
    
    if (nkmers == 0) {
      (*bad_reads)++;
    }
    else {
      Element * current_node  = NULL;

      
      for(i=0;i<windows->nwindows;i++){ //for each window
	KmerSlidingWindow * current_window = &(windows->window[i]);
	
	for(j=0;j<current_window->nkmers;j++){ //for each kmer in window

	  	  current_node = hash_table_find(element_get_key(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph);	  	 

	  if (current_node){
	    db_node_set_status(current_node, exists_in_reference);
	  }
	  
	}
	
      }
    }
   
    if (full_entry == false){
      shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
    }

    prev_full_entry = full_entry;

  }
  
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);

}





// This is just like get_sliding_windows_from_sequence is seq.c, but this one breaks the window also when a kmer is not in the graph.
//The first argument - seq - is a C string in A,C,G,T format
//The second argument - quality - is a string of qualities for the sequence, one byte per base.
//quality cutoff argument defines the threshold for quality
//return total number of kmers read
//The third argument - length - is the length in bases of the sequence.
//return total number of kmers read
int get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(char * seq,  char * qualities, int length, char quality_cut_off, 
										  KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph)  
{


  short kmer_size = db_graph->kmer_size;

  char first_kmer[kmer_size+1];
  first_kmer[kmer_size]='\0';

  BinaryKmer tmp_bin_kmer;
  BinaryKmer tmp_bin_kmer2;
      


  int i=0; //current index
  int count_kmers = 0;

  if (seq == NULL){
    fputs("in get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph, seq is NULL\n",stderr);    
    exit(1);
  }

  if (length < kmer_size || max_windows == 0 || max_kmers == 0){
    return 0;
  }

  int index_windows = 0;
  
  //loop over the bases in the sequence
  //index i is the current position in input sequence -- it nevers decreases. 
  
  do{

    //build first kmer, ie a stretch of kmer_size good qualities bases
    int j = 0; //count how many good bases

    while ((i<length) && (j<kmer_size)){

      //collects the bases in the first kmer
      first_kmer[j] = seq[i];

      if ((char_to_binary_nucleotide(seq[i]) == Undefined) || 
	  (quality_cut_off!=0 && qualities[i]<= quality_cut_off)){
	j=0; //restart the first kmer 
      }
      else{
	j++;
      }

      i++; 
    }    

    if (j==kmer_size){ //ie we did not parse the entire sequence looking for a single good kmer, the first kmer
      
      if (  hash_table_find(element_get_key(seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer), kmer_size, &tmp_bin_kmer2), db_graph) == NULL )
	{
	  j=0;
	  i=i-(kmer_size-1); // first kmer may be bad because of the first base. Start again from just after that
	  continue; //want to restart building a first kmer
	}

      count_kmers++;

      //new sliding window
      if (index_windows>=max_windows){
	  fputs("number of windows is bigger than max_windows",stderr);
	  exit(1);
	}

      KmerSlidingWindow * current_window =&(windows->window[index_windows]);

      int index_kmers = 0;
      //do first kmer
      seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer);
      binary_kmer_assignment_operator(current_window->kmer[index_kmers] , tmp_bin_kmer);

      //do the rest --
      index_kmers++;
    
      while(i<length){
	
	if (index_kmers>=max_kmers){
	  fputs("number of kmers is bigger than max_kmers\n",stderr);
	  exit(1);
	}

	Nucleotide current_base = char_to_binary_nucleotide(seq[i]);
	BinaryKmer tmp_next_kmer;

	//set to previous kmer, then shift an add new base
	binary_kmer_assignment_operator(tmp_next_kmer, current_window->kmer[index_kmers-1]);
	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&tmp_next_kmer, current_base, kmer_size);

	if (  (current_base == Undefined) 
	      ||
	      (quality_cut_off!=0 && qualities[i]<= quality_cut_off)
	      )
	  {
	    break;
	  }
	else if (hash_table_find(element_get_key(&tmp_next_kmer, kmer_size, &tmp_bin_kmer2), db_graph) == NULL)
	  {
	    i++;
	    break;
	  } 
	binary_kmer_assignment_operator(current_window->kmer[index_kmers], tmp_next_kmer);
	index_kmers++;
	count_kmers++;
	i++;
      }

      current_window->nkmers = index_kmers; 
      index_windows++;
            
    }
  } while (i<length);
 

  windows->nwindows = index_windows;

  return count_kmers;



}


// This is just like get_sliding_windows_from_sequence is seq.c, but this one breaks a window if any of the kmers is not in the graph
// OR if any of the edges is not in the graph. ie entire window and edges must lie in graph

//The first argument - seq - is a C string in A,C,G,T format
//The second argument - quality - is a string of qualities for the sequence, one byte per base.
//quality cutoff argument defines the threshold for quality
//return total number of kmers read
//The third argument - length - is the length in bases of the sequence.
//return total number of kmers read
int get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(char * seq,  char * qualities, int length, char quality_cut_off, 
										     KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph)  
{
  short kmer_size = db_graph->kmer_size;

  char first_kmer[kmer_size+1];
  first_kmer[kmer_size]='\0';

  BinaryKmer tmp_bin_kmer;
  BinaryKmer tmp_bin_kmer2;
  BinaryKmer tmp_bin_kmer3;
      


  int i=0; //current index
  int count_kmers = 0;

  if (seq == NULL){
    fputs("in get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph, seq is NULL\n",stderr);    
    exit(1);
  }

  if (length < kmer_size || max_windows == 0 || max_kmers == 0){
    return 0;
  }

  int index_windows = 0;
  
  //loop over the bases in the sequence
  //index i is the current position in input sequence -- it nevers decreases. 
  
  do{ //first do

    //build first kmer, ie a stretch of kmer_size good qualities bases
    int j = 0; //count how many good bases

    while ((i<length) && (j<kmer_size)){

      //collects the bases in the first kmer
      first_kmer[j] = seq[i];

      if ((char_to_binary_nucleotide(seq[i]) == Undefined) || 
	  (quality_cut_off!=0 && qualities[i]<= quality_cut_off)){
	j=0; //restart the first kmer 
      }
      else{
	j++;
      }

      i++; 
    }    

    if (j==kmer_size){ //ie we did not parse the entire sequence looking for a single good kmer, the first kmer
      
      if (  hash_table_find(element_get_key(seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer), kmer_size, &tmp_bin_kmer2), db_graph) == NULL )
	{
	  j=0;
	  i=i-(kmer_size-1); // first kmer may be bad because of the first base. Start again from just after that
	  continue; //want to restart building a first kmer       

	  // give up
	  //i=length;
	  //index_windows=0;
	  //break;
	  
	}

      count_kmers++;

      //new sliding window
      if (index_windows>=max_windows){
	  fputs("number of windows is bigger than max_windows",stderr);
	  exit(1);
	}

      KmerSlidingWindow * current_window =&(windows->window[index_windows]);

      int index_kmers = 0;
      //do first kmer
      seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer);
      binary_kmer_assignment_operator(current_window->kmer[index_kmers] , tmp_bin_kmer);

      //do the rest --
      index_kmers++;
    
      while(i<length){
	
	if (index_kmers>=max_kmers){
	  fputs("number of kmers is bigger than max_kmers\n",stderr);
	  exit(1);
	}

	Nucleotide current_base = char_to_binary_nucleotide(seq[i]);
	BinaryKmer tmp_curr_kmer;

	//set to previous kmer, then shift an add new base
	binary_kmer_assignment_operator(tmp_curr_kmer, current_window->kmer[index_kmers-1]);
	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&tmp_curr_kmer, current_base, kmer_size);

	dBNode* test_prev_node = hash_table_find(element_get_key(&(current_window->kmer[index_kmers-1]) , kmer_size, &tmp_bin_kmer2), db_graph);
	dBNode* test_curr_node = hash_table_find(element_get_key(&tmp_curr_kmer,                          kmer_size, &tmp_bin_kmer3), db_graph);


	if (  (current_base == Undefined) 
	      ||
	      (quality_cut_off!=0 && qualities[i]<= quality_cut_off)
	      )
	  {
	    //give up
	    //i=length;
	    //index_windows=0;
	    break;
	  }
	else if  (test_curr_node == NULL)
	  {
	    //give up
	    //i=length;
	    //index_windows=0;
	    i++;
	    break;
	  } 
	else //check edge exists between prev kmer and this one   
	  {
	    Orientation prev_orientation = db_node_get_orientation(&(current_window->kmer[index_kmers-1]),test_prev_node, kmer_size);
	    
	    if (! (db_node_edge_exist(test_prev_node, current_base, prev_orientation) ) )
	      {
		//give up
		//i=length;
		i=i-(kmer_size-1);
		//i++;
		//index_windows=0;
		break;
	      }

	  }
	binary_kmer_assignment_operator(current_window->kmer[index_kmers], tmp_curr_kmer);
	index_kmers++;
	count_kmers++;
	i++;
      }

      current_window->nkmers = index_kmers; 
      index_windows++;
            
    }
  } while (i<length);
 

  windows->nwindows = index_windows;

  return count_kmers;



}




void read_fastq_and_print_reads_that_lie_in_graph(FILE* fp, FILE* fout, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
						  long long * bad_reads, int max_read_length, dBGraph * db_graph,
						  boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index)
{
    //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  
  int seq_length=0;
  short kmer_size = db_graph->kmer_size;

  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);

  char tmpseq[db_graph->kmer_size+1];
  tmpseq[db_graph->kmer_size]='\0';
  tmpseq[0]='\0';
  

  BinaryKmer tmp_kmer;
  boolean full_entry = true;
  boolean prev_full_entry = true;

  int entry_length;

  while ((entry_length = file_reader(fp,seq,max_read_length, full_entry, &full_entry))){

    if (DEBUG){
      printf ("\nsequence %s\n",seq->seq);
    }
    
    int i,j;
    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));

    
    //use quality cutoff of 0, arg 4 below
    // int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,0,db_graph->kmer_size,windows,max_windows, max_kmers);
    int nkmers = get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq,seq->qual,entry_length,0,
											       windows,max_windows, max_kmers, db_graph);
    
    if (nkmers == 0) 
      {
	(*bad_reads)++;
      }
    else 
      {
	Element * current_node  = NULL;
	//we will throw out any window that does not have ALL kmers in the graph
	for(i=0;i<windows->nwindows;i++)
	  { //for each window
	    KmerSlidingWindow * current_window = &(windows->window[i]);
	    
	    boolean all_kmers_in_this_window_are_in_graph=true;
	    

	    for(j=0;j<current_window->nkmers;j++)
	      { //for each kmer in window
		
		current_node = hash_table_find(element_get_key(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph);	  	 
		if (current_node==NULL)
		  {
		    all_kmers_in_this_window_are_in_graph=false;
		  }
	      }

	    
	    if (all_kmers_in_this_window_are_in_graph==true)
	      {
		//print out this window as a "read". If this read has many windows, we will print each as a separate read (provided they lie in the graph)
		if (is_for_testing == false)
		  {
		    fprintf(fout, "> %s part %d \n", seq->name, i );
		    fprintf(fout, "%s", binary_kmer_to_seq(&(current_window->kmer[0]), db_graph->kmer_size, tmpseq) );
		    for(j=1;j<current_window->nkmers;j++){ 
		      fprintf(fout, "%c", binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&(current_window->kmer[j]))) );
		    }
		    fprintf(fout, "\n");
		  }
		else
		  {
		    for_test_array_of_clean_reads[*for_test_index][0]='\0';
		    strcat(for_test_array_of_clean_reads[*for_test_index], binary_kmer_to_seq(&(current_window->kmer[0]), db_graph->kmer_size, tmpseq) );
		    for(j=1;j<current_window->nkmers;j++){ 
		      char tmp_ch[2];
		      tmp_ch[0]=binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&(current_window->kmer[j])));
		      tmp_ch[1]='\0';
		      strcat(for_test_array_of_clean_reads[*for_test_index], tmp_ch );
		    }
		    *for_test_index=*for_test_index+1;
		  }
	      }
	  //else
	  //  {
	  //  }
	  
	  }
      }
    
    if (full_entry == false){
      shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
    }

    prev_full_entry = full_entry;

  }
  
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);

}


void read_fastq_and_print_subreads_that_lie_in_graph_breaking_at_edges_or_kmers_not_in_graph(FILE* fp, FILE* fout,
											     int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, 
														 boolean * full_entry), 
											     long long * bad_reads, int max_read_length, dBGraph * db_graph,
											     boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index)
{
    //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  
  int seq_length=0;
  short kmer_size = db_graph->kmer_size;

  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);

  char tmpseq[db_graph->kmer_size+1];
  tmpseq[db_graph->kmer_size]='\0';
  tmpseq[0]='\0';
  

  //BinaryKmer tmp_kmer;
  boolean full_entry = true;
  boolean prev_full_entry = true;

  int entry_length;

  while ((entry_length = file_reader(fp,seq,max_read_length, full_entry, &full_entry))){

    if (DEBUG){
      printf ("\nsequence %s\n",seq->seq);
    }
    
    int i,j;
    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));

    
    //use quality cutoff of 0, arg 4 below
    // int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,0,db_graph->kmer_size,windows,max_windows, max_kmers);
    //    int nkmers = get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq,seq->qual,entry_length,0,
    //											       windows,max_windows, max_kmers, db_graph);
    
    int nkmers = get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, entry_length, 0,
												  windows, max_windows, max_kmers, db_graph);
    

    if (nkmers == 0) 
      {
	(*bad_reads)++;
      }
    else 
      {
	//Element * current_node  = NULL;
	
	for(i=0;i<windows->nwindows;i++)
	  { //for each window
	    KmerSlidingWindow * current_window = &(windows->window[i]);
	    
	    //print out this window as a "read". If this read has many windows, we will print each as a separate read (provided they lie in the graph)
	    if (is_for_testing == false)
	      {
		fprintf(fout, "> %s part %d \n", seq->name, i );
		fprintf(fout, "%s", binary_kmer_to_seq(&(current_window->kmer[0]), db_graph->kmer_size, tmpseq) );
		for(j=1;j<current_window->nkmers;j++){ 
		  fprintf(fout, "%c", binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&(current_window->kmer[j]))) );
		}
		fprintf(fout, "\n");
	      }
	    else
	      {
		for_test_array_of_clean_reads[*for_test_index][0]='\0';
		strcat(for_test_array_of_clean_reads[*for_test_index], binary_kmer_to_seq(&(current_window->kmer[0]), db_graph->kmer_size, tmpseq) );
		for(j=1;j<current_window->nkmers;j++){ 
		  char tmp_ch[2];
		  tmp_ch[0]=binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&(current_window->kmer[j])));
		  tmp_ch[1]='\0';
		  strcat(for_test_array_of_clean_reads[*for_test_index], tmp_ch );
		}
		*for_test_index=*for_test_index+1;
	      }
	    
		    
	  }
      }
	
	
    if (full_entry == false){
      shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
    }
    
    prev_full_entry = full_entry;
    
  }
  
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
  
}



void read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(char* f_name, dBGraph* db_graph)
{

  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry)
  {
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }

  int max_read_length = 3000;


  FILE* fptr = fopen(f_name, "r");
  if (fptr == NULL)
    {
      printf("cannot open chromosome fasta file:%s\n",f_name);
      exit(1);
    }

  long long bad_reads;
  read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(fptr, &file_reader, &bad_reads,  0 , max_read_length, db_graph);


}


void read_all_ref_chromosomes_and_mark_graph(dBGraph* db_graph)
{

  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.1.fa.short_reads" , db_graph);
  printf("Loaded chromosome 1\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.2.fa.short_reads" , db_graph);
  printf("Loaded chromosome 2\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.3.fa.short_reads" , db_graph);
  printf("Loaded chromosome 3\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.4.fa.short_reads" , db_graph);
  printf("Loaded chromosome 4\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.5.fa.short_reads" , db_graph);
  printf("Loaded chromosome 5\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.6.fa.short_reads" , db_graph);
  printf("Loaded chromosome 6\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.7.fa.short_reads" , db_graph);
  printf("Loaded chromosome 7\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.8.fa.short_reads" , db_graph);
  printf("Loaded chromosome 8\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.9.fa.short_reads" , db_graph);
  printf("Loaded chromosome 9\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.10.fa.short_reads" , db_graph);
  printf("Loaded chromosome 10\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.11.fa.short_reads" , db_graph);
  printf("Loaded chromosome 11\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.12.fa.short_reads" , db_graph);
  printf("Loaded chromosome 12\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.13.fa.short_reads" , db_graph);
  printf("Loaded chromosome 13\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.14.fa.short_reads" , db_graph);
  printf("Loaded chromosome 14\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.15.fa.short_reads" , db_graph);
  printf("Loaded chromosome 15\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.16.fa.short_reads" , db_graph);
  printf("Loaded chromosome 16\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.17.fa.short_reads" , db_graph);
  printf("Loaded chromosome 17\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.18.fa.short_reads" , db_graph);
  printf("Loaded chromosome 18\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.19.fa.short_reads" , db_graph);
  printf("Loaded chromosome 19\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.20.fa.short_reads" , db_graph);
  printf("Loaded chromosome 20\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.21.fa.short_reads" , db_graph);
  printf("Loaded chromosome 21\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.22.fa.short_reads" , db_graph);
  printf("Loaded chromosome 22\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.X.fa.short_reads" , db_graph);
  printf("Loaded chromosome X\n");
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/split/clean/Homo_sapiens.NCBI36.52.dna.chromosome.Y.fa.short_reads" , db_graph);
  printf("Loaded chromosome Y\n");
  
}



