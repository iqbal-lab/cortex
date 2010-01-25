/*
  routines to load files into dB Graph that is aware of multiple people
 */

#include <binary_kmer.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <dB_graph.h>
#include <seq.h>
#include <file_reader.h>
#include <global.h>
#include <string.h>
#include <dB_graph_supernode.h>

int MAX_FILENAME_LENGTH=500;
int MAX_READ_LENGTH=10000;

long long load_seq_data_into_graph_of_specific_person_or_pop(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
							     long long * bad_reads, char qualiy_cut_off, long long* dup_reads,
							     int max_read_length, boolean remove_dups_single_endedly, boolean break_homopolymers, int homopolymer_cutoff, 
							     dBGraph* db_graph, EdgeArrayType type, int index);

long long load_paired_end_seq_into_graph_of_specific_person_or_pop(FILE* fp1, FILE* fp2, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry),
								   long long * bad_reads, char quality_cut_off, int max_read_length, long long* dup_reads, boolean remove_dups, 
								   boolean break_homopolymers, int homopolymer_cutoff, dBGraph * db_graph, EdgeArrayType type, int index);



long long load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * bad_reads,  char quality_cut_off, long long* dup_reads, int max_read_length, 
									     boolean remove_duplicates_single_endedly, boolean break_homopolymers, int homopolymer_cutoff,
									     dBGraph* db_graph, EdgeArrayType type, int index)
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

  ret =  load_seq_data_into_graph_of_specific_person_or_pop(fp,&file_reader,bad_reads,quality_cut_off,dup_reads, max_read_length,remove_duplicates_single_endedly, 
							    break_homopolymers, homopolymer_cutoff, db_graph, type, index);

  fclose(fp);

  return ret;
}

long long load_paired_fastq_from_filenames_into_graph_of_specific_person_or_pop(char* filename1, char* filename2, 
										long long * bad_reads,  char quality_cut_off, int max_read_length, 
										long long* dup_reads, boolean remove_duplicates, boolean break_homopolymers, int homopolymer_cutoff, 
										dBGraph* db_graph, EdgeArrayType type, int index )
{

  printf("Start loading %s and %s\n", filename1, filename2);
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

  long long ret =  load_paired_end_seq_into_graph_of_specific_person_or_pop(fp1, fp2, &file_reader,bad_reads,quality_cut_off,max_read_length, dup_reads, remove_duplicates, 
									    break_homopolymers, homopolymer_cutoff, db_graph, type, index);
  fclose(fp1);
  fclose(fp2);

  printf("Finished loading. Total duplicate reads in %s and %s: %lld\n", filename1, filename2, *dup_reads);
  return ret;
}


//this routine supports big fasta entries (chromosome length for example)
long long load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * bad_reads, long long* dup_reads, int max_chunk_length, 
									     boolean remove_duplicates_single_endedly, boolean break_homopolymers, int homopolymer_cutoff, 
									     dBGraph* db_graph, EdgeArrayType type, int index)
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
  ret =  load_seq_data_into_graph_of_specific_person_or_pop(fp,&file_reader,bad_reads,0,dup_reads, max_chunk_length,remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
								db_graph, type, index);
  
  fclose(fp);


  return ret;
}



//do not export the folloiwing internal function
void load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(KmerSlidingWindowSet * windows, boolean* prev_full_ent, 
											     boolean* full_ent, long long* seq_len, boolean mark_read_starts, 
											     dBGraph* db_graph, EdgeArrayType type, int index)
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
	    db_node_update_coverage(current_node,type, index, 1);
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
		   binary_kmer_to_seq(binary_kmer_reverse_complement(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph->kmer_size,kmer_seq2),
		   db_node_get_coverage(current_node, type, index));
	  }
	  
	  if (j>0){
	    
	    if (previous_node == NULL){
	      printf("PROBLEM i:%i j:%i seq_length:%qd nkmers:%i prev_full_entry:%s\n",i,j,*(seq_len),current_window->nkmers,*prev_full_ent == true ? "true" : "false");
	      fputs("file_reader: problem - prev kmer not found\n",stderr);
	      exit(1);
	    }
	    else{ //this is the point at which the new element/node gets associated with the specific person
	      db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size, type, index);
	    }
	  }
	  previous_node = current_node;
	  previous_orientation = current_orientation;
	  
	}
      }
     
}


//pass in a single kmer sliding window and the Sequence* it was derived from. Will find the nodes correspinding to this seqeunce
//and put them in array. Also will check that edges exist as expected from the Sequence*
void load_kmers_from_sliding_window_into_array(KmerSlidingWindow* kmer_window, Sequence* seq, dBGraph* db_graph, dBNode** array_nodes, Orientation* array_orientations, 
					       int max_array_size, 
					       boolean require_nodes_to_lie_in_given_colour, int colour)

{

      Element * current_node  = NULL;
      Element * previous_node  = NULL;
      Orientation current_orientation=forward;
      Orientation previous_orientation=forward;
      Nucleotide current_base;
      BinaryKmer tmp_kmer;
      int j;
	
      if (kmer_window->nkmers>max_array_size)
	{
	  printf("Cannot load_kmers_from_sliding_window_into_array as max_array_size %d < number f kmers in the window %d\n", max_array_size, kmer_window->nkmers);
	  exit(1);
	}

      for(j=0;j<kmer_window->nkmers;j++){ //for each kmer in window
	current_node = hash_table_find(element_get_key(&(kmer_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph);	  
	if (current_node == NULL){
	  fputs("load_kmers_from_sliding_window_into_array: problem - current kmer not found\n",stderr);
	  exit(1);
	}
	if ( (require_nodes_to_lie_in_given_colour==true) && 
	     (db_node_is_this_node_in_this_person_or_populations_graph(current_node, individual_edge_array, colour)==false) )
	  {
	    printf("This current node does not exist in colour %d\n", colour);
	    exit(1);
	  }

	  
	//if not adding the checks for edges, then can dlete all the orientation stuff.
	current_orientation = db_node_get_orientation(&(kmer_window->kmer[j]),current_node, db_graph->kmer_size);

	//add to array
	array_nodes[j]        = current_node;
	array_orientations[j] = current_orientation;

	if (j>0)
	  {
	    current_base = char_to_binary_nucleotide(seq->seq[db_graph->kmer_size -1+j ]);
	    if (previous_node == NULL)
	      {
		printf("PROBLEM j:%i nkmers:%i \n",j, kmer_window->nkmers);
		printf("file_reader: we have a problem - prev kmer not found\n");
		exit(1);
	      }
	    else if (! (db_node_edge_exist(previous_node, current_base, previous_orientation, individual_edge_array, colour) ) )
	      { 
		printf("Missing edge");
		exit(1);
	      }
	  }
	previous_node = current_node;
	previous_orientation = current_orientation;
	  
      }
}





//returns length of sequence loaded
//penultimate argument is to specify whether to discard potential PCR duplicate reads single-endedly - ie if read starts at same kmer
// as a previous read, then discard it. This is a pretty harsh filter, and if ppssible, prefer to use paired end info.
// So in general when calling this function, would expect that boolean remove_dups_single_endedly to be set to false, unless you know you have low coverage, so have
// low probability of two reads starting at the same point.
long long load_seq_data_into_graph_of_specific_person_or_pop(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
							     long long * bad_reads, char quality_cut_off, long long* dup_reads, 
							     int max_read_length, boolean remove_dups_single_endedly, 
							     boolean break_homopolymers, int homopolymer_cutoff, dBGraph * db_graph, EdgeArrayType type, int index){

 
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
    

    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));
   
    //printf("Length %qd %i %s %i %i\n",seq_length,entry_length, seq->name, seq->start, seq->end);

 
    int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,
						   entry_length,quality_cut_off,db_graph->kmer_size,windows,max_windows, 
						   max_kmers, break_homopolymers, homopolymer_cutoff);

    //printf("number kmers:%i %i\n",nkmers,windows->nwindows);


    
    if (nkmers == 0) {
      (*bad_reads)++;
    }
    else {


      //we are assuming that if you want to remove dups, you deal with one sample at a time, in a single coloured graph.
      //I don't want to comtempate multiple colours but only one type of read_start label. So I'm happy to use hash_table_find here
      //and not check if it is in this person's/colour of the graph
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
		  //printf("Discard duplicate single ended read %s\n", seq->name);
		  continue;
		}
	    }
	}
      
      //if remove_dups_single_endedly==true/false, then we are passing in true/false for the mark_read_starts argument
      load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows, &prev_full_entry, &full_entry, &seq_length, remove_dups_single_endedly, 
											      db_graph, type, index);

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
void paired_end_sequence_core_loading_loop_of_specific_person_or_pop(FILE* fp1 , FILE* fp2, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry),
								     Sequence* seq1, Sequence* seq2, char quality_cut_off, int max_read_length, int max_kmers, int max_windows, 
								     KmerSlidingWindowSet * windows1, KmerSlidingWindowSet * windows2, long long* seq_len, long long* dup_reads, long long* bad_reads, 
								     boolean remove_dups, boolean break_homopolymers, int homopolymer_cutoff, 
								     dBGraph* db_graph, EdgeArrayType type, int index )
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
						    entry_length1,quality_cut_off,db_graph->kmer_size,windows1,max_windows, max_kmers, break_homopolymers, homopolymer_cutoff);
    int nkmers2 = get_sliding_windows_from_sequence(seq2->seq,seq2->qual,
						    entry_length2,quality_cut_off,db_graph->kmer_size,windows2,max_windows, max_kmers, break_homopolymers, homopolymer_cutoff);


    //check whether first kmer of both reads is a read-start. If yes, then discard read as duplicate.
    BinaryKmer tmp_kmer;

    if ( (remove_dups==true) && (nkmers1!=0) && (nkmers2 !=0) )
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
		//printf("Discard duplicate reads %s and %s\n", seq1->name, seq2->name);
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
	load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows2, &prev_full_entry2, &full_entry2, seq_len, true, db_graph, type, index);
      }
    else if ((nkmers1!=0)&&(nkmers2==0))
      {
	(*bad_reads)++;
	load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows1, &prev_full_entry1, &full_entry1, seq_len, true, db_graph, type, index);
      }
    else 
      {
	load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows1, &prev_full_entry1, &full_entry1, seq_len, true, db_graph, type, index);
	load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows2, &prev_full_entry2, &full_entry2, seq_len, true, db_graph, type, index);
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
  

//do not export
//returns length of sequence loaded
long long load_paired_end_seq_into_graph_of_specific_person_or_pop(FILE* fp1, FILE* fp2, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
								   long long * bad_reads, char quality_cut_off, int max_read_length, long long* dup_reads, boolean remove_dups, 
								   boolean break_homopolymers, int homopolymer_cutoff, dBGraph * db_graph, EdgeArrayType type, int index){

 


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
  paired_end_sequence_core_loading_loop_of_specific_person_or_pop(fp1, fp2, file_reader, seq1, seq2, quality_cut_off,  max_read_length, max_kmers, max_windows, 
								  windows1, windows2,  &(seq_length), dup_reads, bad_reads, remove_dups, break_homopolymers, homopolymer_cutoff,db_graph, type, index);
  


  free_sequence(&seq1);
  free_sequence(&seq2);
  binary_kmer_free_kmers_set(&windows1);
  binary_kmer_free_kmers_set(&windows2);
  return seq_length;    
}



//assume we have two lists of equal length, of mate fastq files in the same order in each file
long long load_list_of_paired_end_fastq_into_graph_of_specific_person_or_pop(char* list_of_left_mates, char* list_of_right_mates, char quality_cut_off, int max_read_length, 
									     long long* bad_reads, long long* num_dups, int* count_file_pairs, boolean remove_dups, 
									     boolean break_homopolymers, int homopolymer_cutoff,dBGraph* db_graph, EdgeArrayType type, int index)
{

  FILE* f1 = fopen(list_of_left_mates, "r");
  if (f1==NULL)
    {
      printf("Cannot open left mate file %s - exit\n", list_of_left_mates);
      exit(1);
    }
  FILE* f2 = fopen(list_of_right_mates, "r");
  if (f2==NULL)
    {
      printf("Cannot open right mate file %s - exit\n", list_of_right_mates);
      exit(1);
    }

  long long seq_length = 0;
  char filename1[MAX_FILENAME_LENGTH+1];
  char filename2[MAX_FILENAME_LENGTH+1];
  filename1[0]='\0';
  filename2[0]='\0';


  while (feof(f1) ==0)
    {
      if (fgets(filename1,MAX_FILENAME_LENGTH, f1) !=NULL)
	{

	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(filename1, '\n')) != NULL)
	    {
	      *p = '\0';
	    }
	  
	  if (fgets(filename2,MAX_FILENAME_LENGTH, f2) ==NULL)
	    {
	      printf("Files %s and %s are not the same length, and they are lists of mate fastqs\n",list_of_left_mates, list_of_right_mates);
	    }
	  if ((p = strchr(filename2, '\n')) != NULL)
	    {
	      *p = '\0';
	    }
	  
	  seq_length += load_paired_fastq_from_filenames_into_graph_of_specific_person_or_pop(filename1, filename2, bad_reads, quality_cut_off, max_read_length, num_dups, remove_dups, 
											      break_homopolymers, homopolymer_cutoff, db_graph, type, index);
	  (*count_file_pairs)++;
	 
	}


    }
    fclose(f1);
    fclose(f2);


    return seq_length;

}








// gets the next number_of_bases_to_load bases from fasta file, and returns them in the array of nodes.
// assumes this file has already been loaded into the graph.
// returns the number of nodes loaded. If this is less than what you asked for, you know it has hit the end of the file.
// We expect this to be used as follows:
// repeated calls of this function load etc into the LAST number_of_bases_to_load places of the relevant arrays

int load_seq_into_array(FILE* chrom_fptr, int number_of_nodes_to_load, int length_of_arrays, 
			dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
			Sequence* seq, KmerSlidingWindow* kmer_window, boolean expecting_new_fasta_entry, dBGraph * db_graph)
{
  
  int offset; //nodes are placed in the array from offset to offset+number_of_nodes_to_load-1, offset is calculated here to ensure they go at the end of the array
              //path_labels and path_string will contain only the edges - ie not include the first kmer.

  int offset_for_filereader=0;//this is the length of seq->seq that is preserved from last time

  if (expecting_new_fasta_entry==false)
    {
      offset = length_of_arrays-number_of_nodes_to_load-1;//the last node of the previous chunk will be found again, and put back in same place in array , hence the -1
      offset_for_filereader=db_graph->kmer_size;

      if (offset<0)
	{
	  printf("Error, offset<0. NOT expecting new fasta entry, length of arrays is %d, num of bases to load is %d, and offset is %d, their difference\n", length_of_arrays, number_of_nodes_to_load, offset);
	  exit(1);
	}
   
    }
  else //expecting a new entry
    {
      offset = length_of_arrays-number_of_nodes_to_load;

      if (offset<0)
        {
          printf("Error, offset<0. YES, expecting new fasta entry, length of arrays is %d, num of bases to load is %d, and offset is %d, defined by length - (bases to load _ kmer_size -1)\n", 
		 length_of_arrays, number_of_nodes_to_load, offset);
	  exit(1);
        }
    }

  long long seq_length=0;
  short kmer_size = db_graph->kmer_size;
  
  boolean full_entry=false;

  int chunk_length;
  int j;

  BinaryKmer marked_kmer; //will have all longlongs in array being ~0
  //  for (j=0; j<NUMBER_OF_BITFIELDS_IN_BINARY_KMER; j++)
  //    {
  //      marked_kmer[j]=~0;
  //    }
  binary_kmer_set_all_bitfields(marked_kmer, ~( (bitfield_of_64bits) 0) );
  


  if (expecting_new_fasta_entry==false)
    {
      //3rd argument is limit set on number of bases in seq before read_sequence_from_fasta returns. We want number_of_nodes_to_load new bases, plus the kmer's woorth of bases already in seq
      chunk_length = read_sequence_from_fasta(chrom_fptr,seq,number_of_nodes_to_load+db_graph->kmer_size, expecting_new_fasta_entry, &full_entry, offset_for_filereader);
    }
  else
    {
      chunk_length = read_sequence_from_fasta(chrom_fptr,seq,number_of_nodes_to_load+db_graph->kmer_size-1, expecting_new_fasta_entry, &full_entry, offset_for_filereader);
    }




  //doesn't matter whether start of fasta entry or not. If YES, then we ignore the first kmer_size bases, as we are only interested in edges between nodes.
  // If NO, then seq has been preloaded with the last k bases from the previous time, which we want to ignore.
  for (j=0; j < number_of_nodes_to_load; j++)
    {
      path_string[offset+j]=seq->seq[db_graph->kmer_size+j];
    }
  path_string[offset+number_of_nodes_to_load]='\0';

    if (DEBUG){
    printf ("\n sequence returned from read_sequence_from_fasta to load_seq_into_array is %s - kmer size: %i - number of bases loaded inc the preassigned ones at start f seq  length: %i \n",seq->seq,db_graph->kmer_size, chunk_length);
    }
  
  j=0;
  seq_length += (long long) chunk_length;
  
  //number of nodes may be less than what we asked for, if we hit the end of the file
  int num_nodes = get_single_kmer_sliding_window_from_sequence(seq->seq,chunk_length,db_graph->kmer_size, kmer_window);
  
  //sanity
  if (expecting_new_fasta_entry==true)
    {      
      if (num_nodes > number_of_nodes_to_load) 
	{
	  printf("Returned more nodes than asked for in load_seq_into_array. Exit/\n");
	  exit(1);
	}
    }
  else
    {
      //we expect number of nodes = num_of_nodes_to_load , except we had to preload seq with kmer_size bases, so that will give us an extra node at the start, which we need to get the first edge 
      if (num_nodes > number_of_nodes_to_load+1) 
	{
          printf("Returned morenodes than asked for inload_seq_into_array. Exit/\n");
          exit(1);

	}
    }
  

  Element * current_node  = NULL;
  Element * previous_node = NULL;
  
  Orientation current_orientation,previous_orientation;
  BinaryKmer tmp_kmer;

  //if num_kmers=0 (ie have hit end of file), then kmer_window->nkmers=0, so will skip this next "for" loop
  for(j=0;j<kmer_window->nkmers;j++)
    { //for each kmer in window

      if ( binary_kmer_comparison_operator(kmer_window->kmer[j], marked_kmer) ) //encoding as  1's in all 64 bits of all bitfields in BineryKmer, 
	                                                                        //a non-kmer - ie anything that would have had an N in it
	{
	  //corresponds to a kmer that contains an N
	  path_nodes[offset+j]        =NULL;
	  path_orientations[offset+j] =forward;

	  //iqbal
	  //if (j>1)
	  //  {  //edges to/from this node don't exist  remember path_nodes[m] is edge from node m to node m+1.
	  //    path_labels[offset+j-1] = Undefined;
	  //    path_labels[offset+j]   = Undefined;
	  //  }
	  previous_node=NULL;
	  continue;
	}
      
      boolean found = false;
      //ZAM DEBUG - next line is surely wrong
      //current_node = hash_table_find_or_insert(element_get_key(&(kmer_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),&found,db_graph);	  
      current_node = hash_table_find(element_get_key(&(kmer_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph);	  

      if (current_node == NULL){
	BinaryKmer tmp_dbg_kmer;
	char tmp_dbg_seq[db_graph->kmer_size];

	printf("Problem in load_seq_into_array - current kmer not found %s\n",
	       binary_kmer_to_seq(element_get_key(&(kmer_window->kmer[j]),db_graph->kmer_size, &tmp_dbg_kmer), db_graph->kmer_size,tmp_dbg_seq));
	exit(1);
      }
      
      current_orientation = db_node_get_orientation(&(kmer_window->kmer[j]),current_node, db_graph->kmer_size);

      path_nodes[offset+j]        = current_node;
      path_orientations[offset+j] = current_orientation;


      if (DEBUG){
	char kmer_seq[db_graph->kmer_size];
	printf("j=%d, Current node kmer is  %s\n",j, binary_kmer_to_seq(&(kmer_window->kmer[j]),db_graph->kmer_size,kmer_seq));
	if (current_orientation==forward)
	  printf("Current orientation is forward\n");
	else
	  printf("Current orientation is reverse");
      }
      
      if (j>0)
	{
	  if (previous_node == NULL)
	    {
	      if (j>0)
		{
		  path_labels[offset+j-1]=Undefined;
		}
	    }
	  else
	    {
	      BinaryKmer previous_k, current_k; 
	      char seq1[kmer_size];
	      char seq2[kmer_size];
	      
	      binary_kmer_assignment_operator(previous_k, previous_node->kmer);
	      binary_kmer_assignment_operator(current_k, current_node->kmer);
	      
	      if (previous_orientation == reverse){
		binary_kmer_assignment_operator(previous_k, *(binary_kmer_reverse_complement(&previous_k,kmer_size, &tmp_kmer)) );
	      }
	      
	      if (current_orientation == reverse){
		binary_kmer_assignment_operator(current_k, *(binary_kmer_reverse_complement(&current_k,kmer_size, &tmp_kmer)) ) ;
	      }
    
	      
	      if (DEBUG){
		printf("Found edge %s -%c-> %s\n",binary_kmer_to_seq(&previous_k,kmer_size,seq1),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&current_k)),binary_kmer_to_seq(&current_k,kmer_size,seq2));
	      }
	      if (j>0)
		{
		  path_labels[offset+j-1]=binary_kmer_get_last_nucleotide(&current_k); //iqbal zam added if j>0 and added a -1
		}
	      
	    }
	}


	  
      previous_node = current_node;
      previous_orientation = current_orientation;
      
    }


  if ( (expecting_new_fasta_entry==true) || (num_nodes==0) )
    {
      return num_nodes;
    }
  else
    {
      return num_nodes-1;// remember the first node is just the same again as the last node of the last batch
    }



}





//takes a filename 
// this file contains a list of filenames, each of these represents an individual (and contains a list of fasta for that individual).
long long load_population_as_fasta(char* filename, long long* bad_reads, dBGraph* db_graph)
{

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("load_population_as_fasta cannot open file:%s\n",filename);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
  }


  char line[MAX_FILENAME_LENGTH+1];

  int total_seq_loaded=0;
  int people_so_far=0;

  while(fgets(line,MAX_FILENAME_LENGTH, fp) !=NULL)
    {

      //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';


      people_so_far++;
      if (people_so_far>NUMBER_OF_INDIVIDUALS_PER_POPULATION)
      {
        printf("This filelist contains too many people for a single population, %d", people_so_far);
	exit(1);
      }

      //printf("About to try and load fasta for this person %s\n",line);

      total_seq_loaded = total_seq_loaded + load_all_fasta_for_given_person_given_filename_of_file_listing_their_fasta_files(line, bad_reads, db_graph, people_so_far-1);


    }

  fclose(fp);
  //printf("Finished loading population, witht total seq loaded %d\n",total_seq_loaded); 
  return total_seq_loaded;


}


//index tells you which person within a population it is
int load_all_fasta_for_given_person_given_filename_of_file_listing_their_fasta_files(char* f_name, long long* bad_reads, dBGraph* db_graph, int index)
{
  FILE* fptr = fopen(f_name, "r");
  if (fptr == NULL)
    {
    printf("cannot open person-specific fasta file:%s\n",f_name);
    exit(1); 
    }

  //file contains a list of fasta file names
  char line[MAX_FILENAME_LENGTH+1];
  
  int total_seq_loaded=0;
  
  boolean remove_duplicates_single_endedly=false;
  long long dup_reads = 0;
  boolean break_homopolymers = false;
  int homopolymer_cutoff=0;
  
  while(fgets(line,MAX_FILENAME_LENGTH, fptr) !=NULL)
    {

      //remove newline from endof line- replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';
      
      total_seq_loaded = total_seq_loaded + 
	load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(line, bad_reads, &dup_reads, MAX_READ_LENGTH, 
									   remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
									   db_graph, individual_edge_array, index);

    }

  fclose(fptr);
  return total_seq_loaded;

}




//takes a filename 
// this file contains a list of filenames, each of these represents an individual (and contains a list of fasta for that individual).
long long load_population_as_fastq(char* filename, long long* bad_reads, char quality_cutoff, dBGraph* db_graph)
{

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("load_population_as_fastq cannot open file:%s\n",filename);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
  }

  char line[MAX_FILENAME_LENGTH+1];

  int total_seq_loaded=0;
  int people_so_far=0;

  while(fgets(line,MAX_FILENAME_LENGTH, fp) !=NULL)
    {


      //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';


      people_so_far++;
      if (people_so_far>NUMBER_OF_INDIVIDUALS_PER_POPULATION)
      {
        printf("This filelist contains too many people for a single population, %d", people_so_far);
	exit(1);
      }

      total_seq_loaded = total_seq_loaded + load_all_fastq_for_given_person_given_filename_of_file_listing_their_fastq_files(line, bad_reads, quality_cutoff, db_graph, people_so_far-1);

      printf("Just loaded person number %d, and now have cumulative total of  %d bases with %lld bad_reads so far\n", people_so_far-1, total_seq_loaded, *bad_reads);
    }

  printf("Finished loading population, witht total seq loaded %d\n",total_seq_loaded); 
  fclose(fp);
  return total_seq_loaded;


}


//index tells you which person within a population it is
int load_all_fastq_for_given_person_given_filename_of_file_listing_their_fastq_files(char* f_name, long long* bad_reads, char quality_cutoff,  dBGraph* db_graph, int index)
{
  FILE* fptr = fopen(f_name, "r");
  if (fptr == NULL)
    {
      //printf("cannot open person-specific fasta file:%s\n",f_name);
    exit(1); 
    }

  //file contains a list of fasta file names
  char line[MAX_FILENAME_LENGTH+1];

  long long dup_reads=0;
  int total_seq_loaded=0;
  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  
  while(fgets(line,MAX_FILENAME_LENGTH, fptr) !=NULL)
    {

      //remove newline from endof line- replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';
      
      total_seq_loaded = total_seq_loaded + 
	load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(line, bad_reads, quality_cutoff, &dup_reads,  MAX_READ_LENGTH, 
									   remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
									   db_graph, individual_edge_array, index);

    }

  fclose(fptr);
  return total_seq_loaded;

}





//returns number of kmers loaded
int load_multicolour_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph)
{
  FILE* fp_bin = fopen(filename, "r");
  int seq_length = 0;
  dBNode node_from_file;
  boolean found;
  int count=0;

  if (fp_bin == NULL){
    printf("load_multicolour_binary_data_from_filename_into_graph cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }
  
  //Go through all the entries in the binary file
  while (db_node_read_multicolour_binary(fp_bin,db_graph->kmer_size,&node_from_file)){
    count++;
    
    //if (count % 100000000 == 0 ){
    // printf("loaded %i\n",count);

    //}
   
    dBNode * current_node  = NULL;
    BinaryKmer tmp_kmer;
    current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size, &tmp_kmer),&found,db_graph);
    
    seq_length+=db_graph->kmer_size;
   
    int i;
    for (i=0; i<NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
      {
	add_edges(current_node,individual_edge_array, i, get_edge_copy(node_from_file, individual_edge_array, i));
	db_node_update_coverage(current_node, individual_edge_array, i, db_node_get_coverage(&node_from_file, individual_edge_array,i));
      }

  }
  
  fclose(fp_bin);
  return seq_length;
}





int load_single_colour_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph, EdgeArrayType type, int index)
{


  FILE* fp_bin = fopen(filename, "r");
  int seq_length = 0;
  dBNode node_from_file;
  boolean found;
  int count=0;

  if (fp_bin == NULL){
    printf("load_single_colour_binary_data_from_filename_into_graph cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }
  
  //Go through all the entries in the binary file
  while (db_node_read_single_colour_binary(fp_bin,db_graph->kmer_size,&node_from_file, type, index))
    {
      count++;

      //if (count % 100000000 == 0 ){
      // printf("loaded %i\n",count);
      //}

      dBNode * current_node  = NULL;
      BinaryKmer tmp_kmer;
      current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size, &tmp_kmer),&found,db_graph);
      
      seq_length+=db_graph->kmer_size;//todo - maybe only increment if had to insert, not if was already in graph?
      
      add_edges(current_node,individual_edge_array, index, get_edge_copy(node_from_file, individual_edge_array, index));
      db_node_update_coverage(current_node, individual_edge_array, index, db_node_get_coverage(&node_from_file, individual_edge_array,index) );
    }
  
  fclose(fp_bin);
  return seq_length;

}




long long load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(char* filename,  dBGraph* db_graph, EdgeArrayType type, int index)
{

  FILE* fptr = fopen(filename, "r");
  if (fptr == NULL)
    {
      printf("cannot open %s which is supposed to list all .ctx files for person with index %d \n",filename, index);
      exit(1); 
    }

  //file contains a list of .ctx filenames, as dumped by the graph/ target (NOT sv_trio)
  char line[MAX_FILENAME_LENGTH+1];
  
  int total_seq_loaded=0;
  
  while(fgets(line,MAX_FILENAME_LENGTH, fptr) !=NULL)
    {

      //remove newline from endof line- replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';
      
      total_seq_loaded = total_seq_loaded + 
	load_single_colour_binary_data_from_filename_into_graph(line, db_graph, individual_edge_array, index);

    }

  fclose(fptr);
  return total_seq_loaded;


  
}






//takes a filename 
// this file contains a list of filenames, each of these represents an individual (and contains a list of binaries for that individual).
long long load_population_as_binaries_from_graph(char* filename, dBGraph* db_graph)
{

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("load_population_as_binaries_from_graph cannot open file:%s\n",filename);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
  }

  char line[MAX_FILENAME_LENGTH+1];

  int total_seq_loaded=0;
  int people_so_far=0;

  while(fgets(line,MAX_FILENAME_LENGTH, fp) !=NULL)
    {
      
      //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';


      people_so_far++;
      if (people_so_far>NUMBER_OF_INDIVIDUALS_PER_POPULATION)
      {
        printf("This filelist contains too many people, %d, remember we have set a population limit of %d in variable NUMBER_OF_INDIVIDUALS_PER_POPULATION", 
	       people_so_far,NUMBER_OF_INDIVIDUALS_PER_POPULATION);
	exit(1);
      }

      total_seq_loaded = total_seq_loaded + 
	load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(line, db_graph, individual_edge_array, people_so_far-1);

    }

  printf("Finished loading population, with total seq loaded %d\n",total_seq_loaded); 
  return total_seq_loaded;


}



void read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
									    int max_read_length, dBGraph * db_graph)
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

  //we will be calling functions that apply to fastq as well as fasta. Since we are loading a ref fasta, we do not want to cut homopolymers,
  // and since we have fastq, the qualit-cutfoff os redundant
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  char quality_cut_off=0;
  

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

    
    //last two arguments mean not to break at homopolymers
    int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,quality_cut_off,db_graph->kmer_size,windows,max_windows, 
						   max_kmers,break_homopolymers, homopolymer_cutoff);
    
    if (nkmers == 0) {
      //(*bad_reads)++;
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
//note - this does not take arguments for homopolymer cutting - this in an interna function for aligning fastq to the graph
int get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(char * seq,  char * qualities, int length, char quality_cut_off, 
										     KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph,
										     EdgeArrayType type, int index)  
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
	    
	    if (! (db_node_edge_exist(test_prev_node, current_base, prev_orientation, type, index) ) )
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
											     EdgeArrayType type, int index,
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
												  windows, max_windows, max_kmers, db_graph, type, index);
    

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

  read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(fptr, &file_reader, max_read_length, db_graph);


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




//returns the number of kmers loaded
int align_next_read_to_graph_and_return_node_array(FILE* fp, int max_read_length, dBNode** array_nodes, Orientation* array_orientations, 
						   boolean require_nodes_to_lie_in_given_colour,
						   int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
						   Sequence* seq, KmerSlidingWindow* kmer_window,dBGraph * db_graph, int colour)
{
  
  boolean full_entry = true;
  boolean prev_full_entry = true;

  //get next read as a C string and put it in seq. Entry_length is the length of the read.
  int entry_length = file_reader(fp,seq,max_read_length,full_entry,&full_entry);
  //turn it into a sliding window 
  int nkmers = get_single_kmer_sliding_window_from_sequence(seq->seq,entry_length, db_graph->kmer_size, kmer_window);
  //work through the sliding window and put nodes into the array you pass in.
  load_kmers_from_sliding_window_into_array(kmer_window, seq, db_graph, array_nodes, array_orientations, 
					    max_read_length-db_graph->kmer_size+1, require_nodes_to_lie_in_given_colour, colour);

  return nkmers;
}


// returns 1 unless hits end of file, when returns 0
// Exits if the sequence read does not exist in db_graph (ie any kmer or edge) in the colour specified as last argument

int read_next_variant_from_full_flank_file(FILE* fptr, int max_read_length,
					   dBNode** flank5p,    Orientation* flank5p_or,    int* len_flank5p,
					   dBNode** ref_allele, Orientation* ref_allele_or, int* len_ref_allele, 
					   dBNode** alt_allele, Orientation* alt_allele_or, int* len_alt_allele, 
					   dBNode** flank3p,    Orientation* flank3p_or,    int* len_flank3p,
					   dBGraph* db_graph, int colour)
{
  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;

    if (new_entry!= true){
      printf("new_entry has to be true for read_next_variant_from_full_flank_file\n");
      exit(1);
    }

    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }
  
  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  
  
  //end of intialisation 


  *len_flank5p = align_next_read_to_graph_and_return_node_array(fptr, max_read_length, flank5p, flank5p_or,  true, file_reader,
							       seq, kmer_window, db_graph, colour);

  if (*len_flank5p==0)
    {
      return 0; //end of file (should never have a zero length read as 5prime flank btw)
    }

  // printf("5p flank: %s, length %d\n", seq->seq, *len_flank5p);


  *len_ref_allele = align_next_read_to_graph_and_return_node_array(fptr, max_read_length, ref_allele, ref_allele_or, true, file_reader,
							       seq, kmer_window, db_graph, colour);
  //printf("ref allele flank: %s, length %d\n", seq->seq, *len_ref_allele);
  
  *len_alt_allele = align_next_read_to_graph_and_return_node_array(fptr, max_read_length, alt_allele, alt_allele_or, true, file_reader,
							       seq, kmer_window, db_graph, colour);
  //printf("alt allele: %s, length %d\n", seq->seq, *len_alt_allele );

  *len_flank3p = align_next_read_to_graph_and_return_node_array(fptr, max_read_length, flank3p, flank3p_or, true, file_reader,
								seq, kmer_window, db_graph, colour);
  //printf("3p flank: %s, length %d\n", seq->seq, *len_flank3p);

  if (*len_flank3p==0)
    {
      printf("Malformed full flank format file. Last seq we got back was %s", seq->seq);
      exit(1);
    }

  return 1;

}
					   

