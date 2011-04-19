/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
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
#include <global.h>
#include <string.h>
#include <dB_graph_supernode.h>
#include <dB_graph_population.h>
#include <file_format.h>
#include <unistd.h>

int MAX_FILENAME_LENGTH=500;
int MAX_READ_LENGTH=10000;//should ONLY be used by test code



void load_seq_data_into_graph_of_specific_person_or_pop(FILE* fp, long long* bases_read, long long* bases_pass_filters_and_loaded,long long** readlen_count_array, 
							int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
							long long * bad_reads, char quality_cut_off, long long* dup_reads, 
							int max_read_length, boolean remove_dups_single_endedly, 
							boolean break_homopolymers, int homopolymer_cutoff, dBGraph * db_graph, EdgeArrayType type, int index);

void  load_paired_end_seq_into_graph_of_specific_person_or_pop(FILE* fp1, FILE* fp2, int (* file_reader)(FILE * fp, Sequence * seq, 
													 int max_read_length,boolean new_entry, boolean * full_entry), 
							       long long* bases_read, long long* bases_pass_filters_and_loaded, long long** readlen_count_array,
							       long long * bad_reads, char quality_cut_off, int max_read_length, long long* dup_reads, boolean remove_dups, 
							       boolean break_homopolymers, int homopolymer_cutoff, dBGraph * db_graph, EdgeArrayType type, int index);


//pass in bases_read to track amount of sequence read in, and bases_pass_filters_and_loaded to see how much passed filters and got into the graph
void load_se_and_pe_filelists_into_graph_of_specific_person_or_pop(boolean se, boolean pe, char* se_f, char* pe_f1, char* pe_f2,
								   long long* bases_read, long long* bases_pass_filters_and_loaded,long long** readlen_count_array,
								   int qual_thresh, boolean remv_dups_se, int remv_dups_pe, 
								   boolean break_homopolymers, int homopol_limit, int ascii_fq_offset, FileFormat format, 
								   int max_read_length, int colour, dBGraph* db_graph)
{

  char filename[MAX_FILENAME_LENGTH];
  int num_single_ended_files_loaded   = 0;
  int num_file_pairs_loaded   = 0;


  long long single_seq_bases_read = 0;
  long long single_seq_bases_loaded = 0;
  long long paired_seq_bases_read = 0;
  long long paired_seq_bases_loaded = 0;


  long long bad_se_reads=0;
  long long dup_se_reads=0;

  long long bad_pe_reads=0;
  long long dup_pe_reads=0;

  //Go through all the files, loading data into the graph

  if (se ==true)
    {
      printf("Load single-ended files\n");

      FILE* se_fp = fopen(se_f, "r");
      if (se_fp==NULL)
	{
	  printf("Cannot open this filelist of SE files: %s\n", se_f);
	  exit(1);
	}


      while (!feof(se_fp))
	{
	  int ret_fscan_se = fscanf(se_fp, "%s\n", filename);
	  if (ret_fscan_se==EOF)
	    {
	      printf("Unable to read from %s\n", se_f);
	      exit(1);
	    }

	  num_single_ended_files_loaded++;
	  
	  if (format==FASTQ)
	    {
	      load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(filename, &single_seq_bases_read, &single_seq_bases_loaded,readlen_count_array,
										 &bad_se_reads, qual_thresh, &dup_se_reads, 
										 max_read_length, 
										 remv_dups_se,break_homopolymers, homopol_limit, ascii_fq_offset, 
										 db_graph, individual_edge_array, colour);
	    }
	  else if (format==FASTA)
	    {
	      load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(filename,&single_seq_bases_read, &single_seq_bases_loaded,
										 &bad_se_reads, &dup_se_reads, max_read_length, 
										 remv_dups_se,break_homopolymers, homopol_limit, 
										 db_graph, individual_edge_array, colour);
	    }
	  else
	    {
	      printf("SE data is being passed in, so Format in must be specified as fastq or fasta\n");
	      exit(1);
	    }
	  
	  
	    printf("\nNum SE files loaded:%i\n\tkmers:%qd\n\tCumulative bad reads:%qd\n\tTotal SE sequence parsed:%qd\nTotal SE sequence passed filters and loaded:%qd\n\tDuplicates removed:%qd\n",
	  	 num_single_ended_files_loaded,hash_table_get_unique_kmers(db_graph),bad_se_reads,single_seq_bases_read, single_seq_bases_loaded, dup_se_reads);

	    
	}
      fclose(se_fp);
    }
  if (pe==true)
    {
      printf("Load paired-end files\n");
      load_list_of_paired_end_files_into_graph_of_specific_person_or_pop(pe_f1, pe_f2, format,
									 qual_thresh, max_read_length, 
									 &paired_seq_bases_read, &paired_seq_bases_loaded,readlen_count_array,
									 &bad_pe_reads, &dup_pe_reads, &num_file_pairs_loaded, 
									 remv_dups_pe, break_homopolymers, homopol_limit, ascii_fq_offset,
									 db_graph, individual_edge_array, colour); 
      
      printf("\nNum PE pairs of files loaded:%i\n\tkmers:%qd\n\tCumulative bad reads:%qd\n\tTotal PE sequence parsed:%qd\nTotal PE sequence passed filters and loaded:%qd\n\tDuplicates removed:%qd\n\n",
      	     num_file_pairs_loaded,hash_table_get_unique_kmers(db_graph),bad_pe_reads,paired_seq_bases_read,paired_seq_bases_loaded, dup_pe_reads);
	  
    }
  


  *bases_read                    = *bases_read                    + single_seq_bases_read   +paired_seq_bases_read;
  *bases_pass_filters_and_loaded = *bases_pass_filters_and_loaded + single_seq_bases_loaded + paired_seq_bases_loaded;
  
  
  hash_table_print_stats(db_graph);
      

  
}


void load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long* bases_read, long long* bases_pass_filters_and_loaded, long long** readlen_count_array,
									long long * bad_reads,  char quality_cut_off, long long* dup_reads, int max_read_length, 
									boolean remove_duplicates_single_endedly, boolean break_homopolymers, int homopolymer_cutoff,
									int fastq_ascii_offset,
									dBGraph* db_graph, EdgeArrayType type, int index)
{
  
  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      puts("new_entry has to be true for fastq\n");
      exit(1);
    }

    return read_sequence_from_fastq(fp,seq,max_read_length,fastq_ascii_offset);
  }

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    fprintf(stderr,"cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }


  load_seq_data_into_graph_of_specific_person_or_pop(fp,bases_read, bases_pass_filters_and_loaded,readlen_count_array,
						     &file_reader,bad_reads,quality_cut_off,dup_reads, max_read_length,remove_duplicates_single_endedly, 
						     break_homopolymers, homopolymer_cutoff, db_graph, type, index);
  
  fclose(fp);

}
 

void  load_paired_end_data_from_filenames_into_graph_of_specific_person_or_pop(char* filename1, char* filename2, FileFormat format,
									       long long* bases_read, long long* bases_pass_filters_and_loaded,long long** readlen_count_array,
									       long long * bad_reads,  char quality_cut_off, int max_read_length, 
									       long long* dup_reads, boolean remove_duplicates, boolean break_homopolymers, int homopolymer_cutoff, 
									       int fastq_ascii_offset,
									       dBGraph* db_graph, EdgeArrayType type, int index )
 {


  printf("Start loading %s and %s\n", filename1, filename2);
  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      printf("When loading paired-end data, we demand that reads not be arbitrarily long. Specifically here, the limit is %d\n", max_read_length);
      exit(1);
    }

    int offset = 0; //because we assume new_entry==true


    if (format==FASTQ)
      {
	return read_sequence_from_fastq(fp,seq,max_read_length,fastq_ascii_offset);
      }
    else if (format==FASTA)
      {
	return read_sequence_from_fasta(fp,seq,max_read_length, new_entry, full_entry, offset);
      }
    else
      {
	printf("Format passed into load_paired_end_data_from_filenames_into_graph_of_specific_person_or_pop must be fasta or fastq");
	exit(1);
      }
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

  load_paired_end_seq_into_graph_of_specific_person_or_pop(fp1, fp2, &file_reader, bases_read, bases_pass_filters_and_loaded, readlen_count_array,
							   bad_reads,quality_cut_off,max_read_length, dup_reads, remove_duplicates, 
							   break_homopolymers, homopolymer_cutoff, db_graph, type, index);
  fclose(fp1);
  fclose(fp2);

  //printf("Finished loading. Total duplicate reads in %s and %s: %lld\n", filename1, filename2, *dup_reads);

}


//this routine supports big fasta entries (chromosome length for example)

void load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long* bases_read, long long* bases_pass_filters_and_loaded,
									long long * bad_reads, long long* dup_reads, int max_chunk_length, 
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

  //pass NULL in instead of an array for read length distribution - do not support getting read-length distribution for fasta
  load_seq_data_into_graph_of_specific_person_or_pop(fp,bases_read, bases_pass_filters_and_loaded,NULL,
							    &file_reader,bad_reads,0,dup_reads, max_chunk_length,remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
							    db_graph, type, index);
  
  fclose(fp);


}



//do not export the folloiwing internal function
// this is called repeatedly in a loop.  read_len_count_array is used to collect statistics on the length of reads,
// and needs to allow for the fact that we CUT reads at N's, low quality bases
void  load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(KmerSlidingWindowSet * windows, boolean* prev_full_ent, 
											      boolean* full_ent, long long* bases_loaded, boolean mark_read_starts, 
											      dBGraph* db_graph, EdgeArrayType type, int index, long long** read_len_count_array)
{
  long long total_bases_loaded=0;

  Element * current_node  = NULL;
  Element * previous_node  = NULL;
  Orientation current_orientation=forward;
  Orientation previous_orientation=forward;
  BinaryKmer tmp_kmer;
  int i,j;
 
  for(i=0;i<windows->nwindows;i++){ //for each window
    KmerSlidingWindow * current_window = &(windows->window[i]);

    //update total bases loaded
    long long length_this_window = (long long) (current_window->nkmers+db_graph->kmer_size-1);
    total_bases_loaded+=length_this_window;

    if (read_len_count_array !=NULL)
      {
	//log that we have loaded a "read" of this length, provided it is not the last window - the last window length we will return, so  it can be cumulated, in case this is a long read
	// this should never happen with fastq,as we have set the file_reader not to allow this
	(*(read_len_count_array[length_this_window]))++;
      }


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
	      printf("PROBLEM i:%i j:%i bases loaded:%qd nkmers:%i prev_full_entry:%s\n",i,j,*(bases_loaded),current_window->nkmers,*prev_full_ent == true ? "true" : "false");
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

  *bases_loaded = *bases_loaded + total_bases_loaded;


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
	  //printf("load_kmers_from_sliding_window_into_array: problem - current kmer not found\n");
	  //exit(1);
	}
	if ( (require_nodes_to_lie_in_given_colour==true) && 
	     (db_node_is_this_node_in_this_person_or_populations_graph(current_node, individual_edge_array, colour)==false) )
	  {
	    printf("This current node does not exist in colour %d\n", colour);
	    exit(1);
	  }

	  
	if (current_node !=NULL)
	  {
	    current_orientation = db_node_get_orientation(&(kmer_window->kmer[j]),current_node, db_graph->kmer_size);
	  }
	else
	  {
	    current_orientation = forward;
	  }
	//add to array
	array_nodes[j]        = current_node;
	array_orientations[j] = current_orientation;
	    
	previous_node = current_node;
	previous_orientation = current_orientation;
	  
      }
}






// remove_dups_single_endedly argument is to specify whether to discard potential PCR duplicate reads single-endedly - ie if read starts at same kmer
// as a previous read, then discard it. This is a pretty harsh filter, and if ppssible, prefer to use paired end info.
// So in general when calling this function, would expect that boolean remove_dups_single_endedly to be set to false, unless you know you have low coverage, so have
// low probability of two reads starting at the same point.
void load_seq_data_into_graph_of_specific_person_or_pop(FILE* fp, long long* bases_read, long long* bases_pass_filters_and_loaded, long long** readlen_count_array,
							int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
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
  
  
  //long long seq_length=0;
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
  long long len_last_window=0; //will repeatedly read a chunk from file, break it into windows. The length of the final window
                             // needs to be added to the length of the first window of the next chunk if this chunk was not the end of an entry
  
  while ((entry_length = file_reader(fp,seq,max_read_length,full_entry,&full_entry))){
   
    //printf("entry %s %i %i %s\n",seq->name,seq->start,seq->end,prev_full_entry?"not connect":"connect");
    
 
   if (DEBUG){
      printf ("\nsequence %s - kmer size: %i - entry length: %i - max kmers:%i\n",seq->seq,db_graph->kmer_size, entry_length, max_kmers);
    }
    

    *bases_read = *bases_read + (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));
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
      load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows, &prev_full_entry, &full_entry, bases_pass_filters_and_loaded, 
											      remove_dups_single_endedly, 
											      db_graph, type, index,readlen_count_array);
      
    }

    if (full_entry == false)
      {
	shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
      }

    prev_full_entry = full_entry;
   

  }
  
 
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);

}



//do not export the following internal function
void paired_end_sequence_core_loading_loop_of_specific_person_or_pop(FILE* fp1 , FILE* fp2, 
								     int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length
											 ,boolean new_entry, boolean * full_entry),
								     Sequence* seq1, Sequence* seq2, char quality_cut_off, int max_read_length, int max_kmers, int max_windows, 
								     KmerSlidingWindowSet * windows1, KmerSlidingWindowSet * windows2, 
								     long long* bases_read, long long* bases_pass_filters_and_loaded, long long** readlen_count_array, 
								     long long* dup_reads, long long* bad_reads, 
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
    
    //bases_read is ampount fo sequence read in from files, not necessarily loaded into graph
    (*bases_read) += (long long) (entry_length1 - (prev_full_entry1==false ? db_graph->kmer_size : 0));
    (*bases_read) += (long long) (entry_length2 - (prev_full_entry2==false ? db_graph->kmer_size : 0));

 
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
		(*dup_reads)++;//increment once for each read in the pair
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
	load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows2, &prev_full_entry2, &full_entry2, 
												bases_pass_filters_and_loaded, true, db_graph, type, index, 
												readlen_count_array);
      }
    else if ((nkmers1!=0)&&(nkmers2==0))
      {
	(*bad_reads)++;
	load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows1, &prev_full_entry1, &full_entry1, 
												bases_pass_filters_and_loaded, true, db_graph, type, index, 
												readlen_count_array);
      }
    else 
      {
	load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows1, &prev_full_entry1, &full_entry1, 
												bases_pass_filters_and_loaded, true, db_graph, type, index, 
												readlen_count_array);
	load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(windows2, &prev_full_entry2, &full_entry2, 
												bases_pass_filters_and_loaded, true, db_graph, type, index, 
												readlen_count_array);
      }

    }

    if (full_entry1 == false){
      printf("We only support paired-end files for fastq, where we mandate that read-length <= max_read length. ie no reference chromosomes\n");
      exit(1);
      //shift_last_kmer_to_start_of_sequence(seq1,entry_length1,db_graph->kmer_size);
    }

    if (full_entry2 == false){
      printf("We only support paired-end files for fastq, where we mandate that read-length <= max_read length. ie no reference chromosomes\n");
      exit(1);
      //      shift_last_kmer_to_start_of_sequence(seq2,entry_length2,db_graph->kmer_size);
    }

    prev_full_entry1 = full_entry1;
    prev_full_entry2 = full_entry2;
   

}
  

//do not export

void  load_paired_end_seq_into_graph_of_specific_person_or_pop(FILE* fp1, FILE* fp2, int (* file_reader)(FILE * fp, Sequence * seq, 
													 int max_read_length,boolean new_entry, boolean * full_entry), 
							       long long* bases_read, long long* bases_pass_filters_and_loaded, long long** readlen_count_array, 
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
								  windows1, windows2,  
								  bases_read, bases_pass_filters_and_loaded, readlen_count_array,
								  dup_reads, bad_reads, remove_dups, break_homopolymers, homopolymer_cutoff,
								  db_graph, type, index);
  


  free_sequence(&seq1);
  free_sequence(&seq2);
  binary_kmer_free_kmers_set(&windows1);
  binary_kmer_free_kmers_set(&windows2);

}



//assume we have two lists of equal length, of mate files in the same order in each file
void load_list_of_paired_end_files_into_graph_of_specific_person_or_pop(char* list_of_left_mates, char* list_of_right_mates, FileFormat format,
									char quality_cut_off, int max_read_length, 
									long long* bases_read, long long* bases_loaded, long long** readlen_count_array, 
									long long* bad_reads, long long* num_dups, int* count_file_pairs, boolean remove_dups, 
									boolean break_homopolymers, int homopolymer_cutoff, int fq_ascii_cutoff, 
									dBGraph* db_graph, EdgeArrayType type, int index)
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
	  
	  load_paired_end_data_from_filenames_into_graph_of_specific_person_or_pop(filename1, filename2, format,
										   bases_read, bases_loaded,readlen_count_array, 
										   bad_reads, quality_cut_off, max_read_length, num_dups, remove_dups, 
										   break_homopolymers, homopolymer_cutoff, fq_ascii_cutoff,
										   db_graph, type, index);
	  (*count_file_pairs)++;
	 
	}


    }
    fclose(f1);
    fclose(f2);


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
  
  Orientation current_orientation=forward;
  Orientation previous_orientation=forward;
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
      
      //boolean found = false;
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





// Only used for test code
//takes a filename 
// this file contains a list of filenames, each of these represents an individual (and contains a list of fasta for that individual).
void load_population_as_fasta(char* filename, long long* bases_read, long long* bases_loaded,  long long* bad_reads, dBGraph* db_graph)
{

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("load_population_as_fasta cannot open file:%s\n",filename);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
  }


  char line[MAX_FILENAME_LENGTH+1];


  int people_so_far=0;

  while(fgets(line,MAX_FILENAME_LENGTH, fp) !=NULL)
    {

      //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';


      people_so_far++;
      if (people_so_far>NUMBER_OF_COLOURS)
      {
        printf("This filelist contains too many people for a single population, %d", people_so_far);
	exit(1);
      }

      //printf("About to try and load fasta for this person %s\n",line);

      load_all_fasta_for_given_person_given_filename_of_file_listing_their_fasta_files(line, bases_read, bases_loaded, bad_reads, db_graph, people_so_far-1);


    }

  fclose(fp);
  //printf("Finished loading population, witht total seq loaded %d\n",total_seq_loaded); 



}


//index tells you which person within a population it is
//bases_read is passed in to find out how much sequence there was in the files read-in.
//bases_loaded is passed in to find out how much sequence passed filters (qual, PCR dup, homopol) and was loaded into the graph
void load_all_fasta_for_given_person_given_filename_of_file_listing_their_fasta_files(char* f_name, long long* bases_read, long long* bases_loaded, 
										      long long* bad_reads, dBGraph* db_graph, int index)
{
  FILE* fptr = fopen(f_name, "r");
  if (fptr == NULL)
    {
    printf("cannot open person-specific fasta file:%s\n",f_name);
    exit(1); 
    }

  //file contains a list of fasta file names
  char line[MAX_FILENAME_LENGTH+1];
  
  
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
      
      load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(line, bases_read, bases_loaded,  bad_reads, &dup_reads, MAX_READ_LENGTH, 
									 remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
									 db_graph, individual_edge_array, index);
      
    }

  fclose(fptr);


}




/*
//takes a filename 
// this file contains a list of filenames, each of these represents an individual (and contains a list of fasta for that individual).
void load_population_as_fastq(char* filename, long long* bases_read, long long bases_loaded, long long* bad_reads, char quality_cutoff, int fastq_ascii_offset, dBGraph* db_graph)
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
      if (people_so_far>NUMBER_OF_COLOURS)
      {
        printf("This filelist contains too many people for a single population, %d", people_so_far);
	exit(1);
      }

      load_all_fastq_for_given_person_given_filename_of_file_listing_their_fastq_files(line, bases_read, bases_loaded, bad_reads, quality_cutoff, 
										       fastq_ascii_offset, db_graph, people_so_far-1);

      printf("Loaded person number %d\nCumulative total bases parsed:%qd\nCumulative total bases passing filters and loaded:%qd\nBad reads:%qd\n", people_so_far-1, *bases_read, *bad_reads);
    }

  //printf("Finished loading population, witht total seq loaded %d\n",total_seq_loaded); 
  fclose(fp);


}
*/






//returns number of kmers loaded*kmer_length
 //array_mean_readlens and array_total_seqs are arrays of length NUMBER_OF_COLOURS, so they can hold the mean read length+total seq in every colour
long long load_multicolour_binary_from_filename_into_graph(char* filename,  dBGraph* db_graph, 
							   int* num_cols_in_loaded_binary, int** array_mean_readlens, long long** array_total_seqs)
{

  //printf("Load this binary - %s\n", filename);
  FILE* fp_bin = fopen(filename, "r");
  long long  seq_length = 0;
  dBNode node_from_file;
  element_initialise_kmer_covgs_edges_and_status_to_zero(&node_from_file);

  boolean found;
  int count=0;

  if (fp_bin == NULL){
    printf("load_multicolour_binary_from_filename_into_graph cannot open file:%s\n",filename);
    exit(1); 
  }


  if (!(check_binary_signature(fp_bin, db_graph->kmer_size, BINVERSION, num_cols_in_loaded_binary, array_mean_readlens, array_total_seqs) ) )
    {
      printf("Cannot load this binary - signature check fails. Wrong max kmer, number of colours, or binary version. Exiting.\n");
      exit(1);
    }


  //always reads the multicol binary into successive colours starting from 0 - assumes the hash table is empty prior to this
  while (db_node_read_multicolour_binary(fp_bin,db_graph->kmer_size,&node_from_file, *num_cols_in_loaded_binary)){
    count++;
    
    dBNode * current_node  = NULL;
    BinaryKmer tmp_kmer;
    //current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size, &tmp_kmer),&found,db_graph);
    current_node = hash_table_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size, &tmp_kmer),db_graph);
    
    seq_length+=db_graph->kmer_size;
   
    int i;
    for (i=0; i<(*num_cols_in_loaded_binary) ; i++)
      {
	add_edges(current_node,individual_edge_array, i, get_edge_copy(node_from_file, individual_edge_array, i));
	db_node_update_coverage(current_node, individual_edge_array, i, db_node_get_coverage(&node_from_file, individual_edge_array,i));
      }

  }
  
  fclose(fp_bin);
  return seq_length;
}



// In some special cases (eg if you have pooled individuals and cleaned the graph, dumped a binary, and then reloaded  that in colour 0,
// and now want to load binaries into colour 1 only if the kmer is in the (cleaned) colour 0),
// then we set only_load_kmers_already_in_hash==true. In this case, we only load edges in that overlap with the cleaned graph, in colour_clean
long long load_single_colour_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph, 
								  int* mean_readlen, long long* total_seq,
								  boolean all_entries_are_unique, EdgeArrayType type, int index,
								  boolean only_load_kmers_already_in_hash, int colour_clean)
{

  if ( (only_load_kmers_already_in_hash==true) && (colour_clean>=NUMBER_OF_COLOURS) )
    {
      printf("Called load_single_colour_binary_data_from_filename_into_graph and specified as clean-colour, colour %d, when this executable is compiled for %d colours only. Exit.\n", colour_clean, NUMBER_OF_COLOURS);
      exit(1);
    }


  //printf("Open single colour binary: %s\n", filename);

  FILE* fp_bin = fopen(filename, "r");
  long long  seq_length = 0;
  dBNode tmp_node;
  element_initialise_kmer_covgs_edges_and_status_to_zero(&tmp_node);

  boolean found;
  int count=0;

  if (fp_bin == NULL){
    printf("load_single_colour_binary_data_from_filename_into_graph cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }



  int num_cols_in_binary;
  if (!(check_binary_signature(fp_bin, db_graph->kmer_size, BINVERSION, &num_cols_in_binary, &mean_readlen, &total_seq) ) )
    {
      printf("Cannot load this binary - fails signature check. Exiting.\n");
      exit(1);
    }

  if (num_cols_in_binary!=1)
    {
      printf("Expecting a single colour binary, but instead this one has %d colours\n. Exiting.\n", num_cols_in_binary);
      exit(1);
    }

  
  //Go through all the entries in the binary file
  // each time you load the info into a temporary node, and load them *** into colour number index ***
  while (db_node_read_single_colour_binary(fp_bin,db_graph->kmer_size,&tmp_node, type, index))
    {
      count++;

      dBNode * current_node  = NULL;
      BinaryKmer tmp_kmer;
      if (only_load_kmers_already_in_hash==false) //normal case
	{
	  if (!all_entries_are_unique)
	    {
	      current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&tmp_node),db_graph->kmer_size, &tmp_kmer),&found,db_graph);
	    }
	  else
	    {
	      current_node = hash_table_insert(element_get_key(element_get_kmer(&tmp_node),db_graph->kmer_size, &tmp_kmer), db_graph);
	    }
	  seq_length+=db_graph->kmer_size;
	  add_edges(current_node,individual_edge_array, index, get_edge_copy(tmp_node, individual_edge_array, index));
	  db_node_update_coverage(current_node, individual_edge_array, index, db_node_get_coverage(&tmp_node, individual_edge_array,index) );
	}
      else
	{//check if node exists in hash already. If yes, then load edge and covg info into the appropriate colour
	  current_node = hash_table_find(element_get_key(element_get_kmer(&tmp_node),db_graph->kmer_size, &tmp_kmer), db_graph);
	  if (current_node !=NULL)
	    {
	      Edges pre_existing_edge = get_edge_copy(*current_node, individual_edge_array, colour_clean);
	      Edges edge_from_binary  = get_edge_copy(tmp_node, individual_edge_array, index);
	      Edges edge_to_load = pre_existing_edge & edge_from_binary;
	      add_edges(current_node,individual_edge_array, index,edge_to_load);
	      db_node_update_coverage(current_node, individual_edge_array, index, db_node_get_coverage(&tmp_node, individual_edge_array,index));
	    }
	}

    }

  fclose(fp_bin);
  return seq_length;

}



//ordinarily, only_load_kmers_already_in_hash==false, and colour_clean is ignored.
// If you have a clean graph in colour 0, and you only want to load nodes from the binaries that overlap with this,
// then set only_load_kmers_already_in_hash==true, and specify colour_clean to be that clean graph colour. Usually this is zero,
// we compile for 2 colours only, and we are loading into colour 1.
long long load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(char* filename,  dBGraph* db_graph, GraphInfo* db_graph_info, 
											   boolean all_entries_are_unique, EdgeArrayType type, int index,
											   boolean only_load_kmers_already_in_hash, int colour_clean)
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
      
      //printf("Load this binary: %s, into this colour : %d\n", line, index);
      int mean_read_len_in_this_binary=0;
      long long total_seq_in_this_binary=0;
      total_seq_loaded += 
	load_single_colour_binary_data_from_filename_into_graph(line, db_graph, &mean_read_len_in_this_binary,&total_seq_in_this_binary,
								all_entries_are_unique, individual_edge_array, index,
								only_load_kmers_already_in_hash, colour_clean);
      all_entries_are_unique=false;

      graph_info_update_mean_readlen_and_total_seq(db_graph_info, index,mean_read_len_in_this_binary, total_seq_in_this_binary);

    }

  fclose(fptr);
  return total_seq_loaded;


  
}






//takes a filename 
// this file contains a list of filenames, each of these represents an individual (and contains a list of binaries for that individual).
// these go into successive colours, starting with first_colour
long long load_population_as_binaries_from_graph(char* filename, int first_colour,boolean about_to_load_first_binary_into_empty_graph, 
						 dBGraph* db_graph, GraphInfo* db_graph_info, boolean only_load_kmers_already_in_hash, int colour_clean)
{

  if ( (about_to_load_first_binary_into_empty_graph==true) && (only_load_kmers_already_in_hash==true) )
    {
      printf("You are trying to load binaries into an empty hash table, but are specifying that they should be compared with the existing hash table. User error\n");
      exit(1);
    }

  //printf("Open this list of colours: %s\n", filename);
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("load_population_as_binaries_from_graph cannot open file:%s\n",filename);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
  }

  char line[MAX_FILENAME_LENGTH+1];

  int total_seq_loaded=0;
  int which_colour=first_colour;



  while(fgets(line,MAX_FILENAME_LENGTH, fp) !=NULL)
    {
      
      //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';


      if (which_colour>NUMBER_OF_COLOURS-1)
      {
        printf("This filelist contains too many people, remember we have set a population limit of %d in variable NUMBER_OF_COLOURS. Cannot load into colour %d", 
	       NUMBER_OF_COLOURS, which_colour);
	exit(1);
      }

      //printf("Open this filelist of binaries, %s,  all corresponding to the same colour:%d\n",
      //	     line, which_colour-1);

      
      total_seq_loaded = total_seq_loaded + 
	load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(line, db_graph,db_graph_info, 
											 about_to_load_first_binary_into_empty_graph, 
											 individual_edge_array, which_colour,
											 only_load_kmers_already_in_hash, colour_clean);
      about_to_load_first_binary_into_empty_graph=false;
      which_colour++;
    }

  //printf("Finished loading population, with total seq loaded %d\n",total_seq_loaded); 
  return total_seq_loaded;


}


// takes a filename (a colour_list)
// this file contains a list of filenames (one per colour), each of these represents an individual (and contains a list of single-colour binaries for that individual).
// Also takes two colour numbers. clean_colour is the clean colour, and in_colour is the colour ALL of these individuals are loaded into, thus:
// First take person 0's list of binaries and load them all into colour in_colour, BUT only load nodes that are already in the hash,
//  and only load those edges that are in the clean colour. Then dump a single-colour binary of colour in_colour,
// with filename = colour name PLUS a suffix added on the end.
void dump_successive_cleaned_binaries(char* filename, int in_colour, int clean_colour, char* suffix, dBGraph* db_graph, GraphInfo* db_graph_info )
{

  if (in_colour==clean_colour)
    {
      printf("In dump_successive_cleaned_binaries You cannot specify the same colour as both clean_colour and in_colour\n");
      exit(1);
    }

  //printf("Open this list of colours: %s\n", filename);
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("dump_successive_cleaned_binaries cannot open file:%s\n",filename);
    exit(1);
  }

  char line[MAX_FILENAME_LENGTH+1];

  int total_seq_loaded=0;

  while(fgets(line,MAX_FILENAME_LENGTH, fp) !=NULL)
    {
      
      //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';

      total_seq_loaded = total_seq_loaded + 
	load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(line, db_graph,db_graph_info, false,
											 individual_edge_array, in_colour,
											 true, clean_colour);
      char outfile[1000];
      outfile[0]='\0';
      sprintf(outfile,"%s_%s.ctx",line,suffix);
      db_graph_dump_single_colour_binary_of_specified_colour(outfile, &db_node_condition_always_true,db_graph,db_graph_info,in_colour);
      //reset that colour:
      db_graph_wipe_colour(in_colour,db_graph);
      graph_info_set_seq(db_graph_info, in_colour,0);
      graph_info_set_mean_readlen(db_graph_info, in_colour, 0);
    }

  fclose(fp);
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




//returns the number of kmers loaded
int align_next_read_to_graph_and_return_node_array(FILE* fp, int max_read_length, dBNode** array_nodes, Orientation* array_orientations, 
						   boolean require_nodes_to_lie_in_given_colour,
						   int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
						   Sequence* seq, KmerSlidingWindow* kmer_window,dBGraph * db_graph, int colour)
{
  
  boolean full_entry = true;

  //get next read as a C string and put it in seq. Entry_length is the length of the read.
  int entry_length = file_reader(fp,seq,max_read_length,full_entry,&full_entry);
  //turn it into a sliding window 
  int nkmers = get_single_kmer_sliding_window_from_sequence(seq->seq,entry_length, db_graph->kmer_size, kmer_window);
  //work through the sliding window and put nodes into the array you pass in. Note this may find NULL nodes if the kmer is not in the graph
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
					   


//array_mean_readlens is an array of length num_cols, giving the mean read length of data loaded into each colour
//array_total_seq is an array of length num_cols, giving the total amount of sequence in each colour (ie sequence loaded, AFTER filtering out by quality, PCR dups, homopolymers etc)
void print_binary_signature(FILE * fp,int kmer_size, int num_cols, int* array_mean_readlens, long long* array_total_seq){
  char magic_number[6];
  int version = BINVERSION;
  
  magic_number[0]='C';
  magic_number[1]='O';
  magic_number[2]='R';
  magic_number[3]='T';
  magic_number[4]='E';
  magic_number[5]='X';
  

  int num_bitfields = NUMBER_OF_BITFIELDS_IN_BINARY_KMER;

  fwrite(magic_number,sizeof(char),6,fp);
  fwrite(&version,sizeof(int),1,fp);
  fwrite(&kmer_size,sizeof(int),1,fp);
  fwrite(&num_bitfields, sizeof(int),1,fp);
  fwrite(&num_cols, sizeof(int), 1, fp);

  int i;
  for (i=0; i<num_cols; i++)
    {
      fwrite(&(array_mean_readlens[i]), sizeof(int), 1, fp);
    }
  for (i=0; i<num_cols; i++)
    {
      fwrite(&(array_total_seq[i]), sizeof(long long), 1, fp);
    }
  fwrite(magic_number,sizeof(char),6,fp);

}


//return yes if signature is consistent
boolean check_binary_signature(FILE * fp,int kmer_size, int bin_version, int* number_of_colours_in_binary, int** array_mean_readlens, long long** array_total_seqs)
{
  int read;
  char magic_number[6];
  boolean ret = false;

  read = fread(magic_number,sizeof(char),6,fp);
  if (read>0 &&
      magic_number[0]=='C' &&
      magic_number[1]=='O' &&
      magic_number[2]=='R' &&
      magic_number[3]=='T' &&
      magic_number[4]=='E' &&
      magic_number[5]=='X' ){

    int version;
    read = fread(&version,sizeof(int),1,fp);
    if (read>0 && version==BINVERSION)
      {

	int kmer_size2;
	read = fread(&kmer_size2,sizeof(int),1,fp);
	if ((read>0) && (kmer_size2 == kmer_size) )
	  {

	    int num_bitfields;
	    read = fread(&num_bitfields,sizeof(int),1,fp);

	    if ( (read>0) && (num_bitfields==NUMBER_OF_BITFIELDS_IN_BINARY_KMER) )
	    
	      {
		int num_cols;
		read = fread(&num_cols,sizeof(int),1,fp);
		
		if ( (read>0) && (num_cols<=NUMBER_OF_COLOURS) && (num_cols>0)  )
		  { 
		    *number_of_colours_in_binary = num_cols;

		    int i;
		    boolean problem=false;
		    for (i=0; (i<num_cols) && (problem==false); i++)
		      {
			int mean_read_len;
			read = fread(&mean_read_len,sizeof(int),1,fp);
			if (read==0)
			  {
			    problem=true;
			  }
			else
			  {
			    *(array_mean_readlens[i]) = mean_read_len;
			  }
		      }
		    if (problem==false)
		      {
			for (i=0; (i<num_cols) && (problem==false); i++)
			  {
			    long long total_seq;
			    read = fread(&total_seq,sizeof(long long),1,fp);
			    if (read==0)
			      {
				problem=true;
			      }
			    else
			      {
				*(array_total_seqs[i]) = total_seq;
			      }
			  }
			
			magic_number[0]='\0';
			magic_number[1]='\0';
			magic_number[2]='\0';
			magic_number[3]='\0';
			magic_number[4]='\0';
			magic_number[5]='\0';
			read = fread(magic_number,sizeof(char),6,fp);
			if ( (problem==false)&&
			     (read>0)&&
			     magic_number[0]=='C' &&
			     magic_number[1]=='O' &&
			     magic_number[2]=='R' &&
			     magic_number[3]=='T' &&
			     magic_number[4]=='E' &&
			     magic_number[5]=='X' )
			  {
			    ret = true;
			  }
			else
			  {
			    printf("Binary is missing read-length/sequence info, or the end-of-header magic number.\n");
			  }
		      }
		    else
		      {
			printf("You are loading  binary with %d colours into a graph with %d colours - incompatible\n",
			       num_cols, NUMBER_OF_COLOURS);
		      }
		  }
		else
		  {
		    printf("Kmer of binary matches the current graph. However this binary was dumped with a different max_kmer size to that of the current graph\n");
		    printf("This binary uses %d bitfields, and current graph uses %d\n", num_bitfields, NUMBER_OF_BITFIELDS_IN_BINARY_KMER);
		  }
	      }
	    else
	      {
		printf("You are loading a binary with kmer=%d into a graph with kmer=%d - incompatible\n", kmer_size2, kmer_size);
	      }
	  }
	else
	  {
	    printf("Binary versions do not match.\n");
	  }
      }
    else
      {
	printf("Binary does not have magic number in header. Corrupt, or not a Cortex binary\n");
      }
    
  }
  else
    {
      printf("Binary fails signature check - was build for different kmer, or with different binary version. (For debug purposes: Magic number is %s and read is %d)\n", magic_number, read);
    }
  return ret;

}

//return true if signature is readable
boolean query_binary(FILE * fp,int* binary_version, int* kmer_size, int* number_of_bitfields, int* number_of_colours_in_binary){
  int read;
  char magic_number[6];
  

  read = fread(magic_number,sizeof(char),6,fp);
  if (read>0 &&
      magic_number[0]=='C' &&
      magic_number[1]=='O' &&
      magic_number[2]=='R' &&
      magic_number[3]=='T' &&
      magic_number[4]=='E' &&
      magic_number[5]=='X' ){

    int version;
    read = fread(&version,sizeof(int),1,fp);
    if (read>0)
      {
	int kmer_size2;
	read = fread(&kmer_size2,sizeof(int),1,fp);
	if (read>0)
	  {
	    int num_bitfields;
	    read = fread(&num_bitfields,sizeof(int),1,fp);

	    if ( read>0 )
	      {
		int num_cols;
		read = fread(&num_cols,sizeof(int),1,fp);
		
		if ( read>0  )
		  { 
		    //int* binary_version, int* kmer_size, int* num_bitfields, int* number_of_colours_in_binar
		    *binary_version = version;
		    *kmer_size = kmer_size2;
		    *number_of_bitfields = num_bitfields;
		    *number_of_colours_in_binary = num_cols;
		  }
		else
		  {
		    return false;
		  }
	      }
	    else
	      {
		return false;
	      }
	  }
	else
	  {
	    return false;
	  }
      }
    else
      {
	return false;
      }
  }
  else
    {
      return false;
    }

  return true;
}



//given a list of filenames, check they all exist, and return the number of them
int get_number_of_files_and_check_existence_from_filelist(char* filelist)
{

  int count=0;

  FILE* fp = fopen(filelist, "r");
  if (fp==NULL)
    {
      printf("Cannot open %s\n", filelist);
      exit(1);
    }

  char file[MAX_FILENAME_LENGTH+1];
  file[0]='\0';
  
  while (feof(fp)==0)
    {
      if (fgets(file, MAX_FILENAME_LENGTH, fp) != NULL)
	{
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(file, '\n')) != NULL)
	    {
	      *p = '\0';
	    }
	  if (access(file, R_OK)==-1)
	    {
	      printf("Cannot access file %s listed in %s\n", file, filelist);
	      exit(1);
	    }
	  else
	    {
	      count++;
	    }

	}
    }
  fclose(fp);

  return count;
}

//assumes we already know how many files in the list, and have preallocated an array to hold the filenames
void get_filenames_from_list(char* filelist, char** array, int len)
{
  int count=0;

  FILE* fp = fopen(filelist, "r");
  if (fp==NULL)
    {
      printf("Cannot open %s\n", filelist);
      exit(1);
    }

  char file[MAX_FILENAME_LENGTH+1];
  file[0]='\0';
  
  while ( (feof(fp)==0) && (count<len) )
    {
      if (fgets(file, MAX_FILENAME_LENGTH, fp) != NULL)
	{
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(file, '\n')) != NULL)
	    {
	      *p = '\0';
	    }

	  strcpy(array[count], file);
	  count++;
	}
    }
  fclose(fp);


  
}

// filename is a list of files, one for each colour. Check they all exists, there are not too many,
// ad that each of them contains a alist of valid binaries.
boolean check_colour_list(char* filename, int kmer)
{
  int num_cols_in_list = get_number_of_files_and_check_existence_from_filelist(filename);


  char** list_cols = malloc( sizeof(char*) * num_cols_in_list);
  if (list_cols==NULL)
    {
      printf("OOM. Give up can't even allocate space for the names of colourss\n");
      exit(1);
    }
  int i;
  for (i=0; i< num_cols_in_list; i++)
    {
      list_cols[i] = malloc(sizeof(char)*500);
      if (list_cols[i]==NULL)
	{
	  printf("OOM. Giveup can't even allocate space for the names of the colour file i = %d\n",i);
	  exit(1);
	}
    }

  get_filenames_from_list(filename, list_cols,num_cols_in_list);

  //first - none of these files should be a cortex binary - typical mistake, and worth checking for
  for (i=0; i< num_cols_in_list; i++)
    {
      FILE* fptr = fopen(list_cols[i], "r");
      if (fptr==NULL)
	{
	  printf("Cannot open %s from list %s\n", list_cols[i], filename);
	  exit(1);
	}
      int bv;
      int k;
      int num_b;
      int num_c;
      int check = query_binary(fptr, &bv, &k, &num_b, &num_c);
      if (check==true)
	{
	  printf("Error with input arguments.\n --colour_list requires a list of files, each of which represent a colour.\n Inside each of those ");
	  printf("should be a list of cortex binaries. \nHowever your colour list %s contains this file %s which is itself a cortex binary not a list of binaries\n", filename, list_cols[i]);
	  exit(1);
	}
      else
	{
	  //check it contains only cortex binaries.

	  int count=0;

	  FILE* fp = fopen(list_cols[i], "r");
	  if (fp==NULL)
	    {
	      printf("Cannot open %s\n", list_cols[i]);
	      exit(1);
	    }

	  char file[MAX_FILENAME_LENGTH+1];
	  file[0]='\0';
	  
	  while ( feof(fp)==0 )
	    {
	      if (fgets(file, MAX_FILENAME_LENGTH, fp) != NULL)
		{
		  //remove newline from end of line - replace with \0
		  char* p;
		  if ((p = strchr(file, '\n')) != NULL)
		    {
		      *p = '\0';
		    }
		  FILE* fp_putative_binary=fopen(file, "r");
		  if (fp_putative_binary==NULL)
		    {
		      printf("Cannot open this file %s\n", file);
		      exit(1);
		    }
		  int bv2;
		  int k2;
		  int num_b2;
		  int num_c2;

		  int check2 = query_binary(fp_putative_binary, &bv2, &k2, &num_b2, &num_c2);
		  fclose(fp_putative_binary);
		  if (check2==false)
		    {
		      printf("Error with input arguments. --colour_list requires a list of files, each of which represent a colour. Inside each of those ");
		      printf("should be a list of cortex binaries. Therefore your colour list %s contains this file %s which should be a list of cortex binaries\n", filename, list_cols[i]);
		      printf("However it contains %s, which is not a cortex binary\n", file);
		      exit(1);
		    }
		  else if (k2 != kmer)
		    {
		      printf("This binary %s is a kmer %d binary, and cannot be loaded into the main graph which has kmer %d\n", file,k2, kmer);
		      exit(1);
		    }
		  else if (num_b2!=NUMBER_OF_BITFIELDS_IN_BINARY_KMER)
		    {
		      printf("Cannot load this binary %s, as it was built with max_kmer_size=%d, whereas the current graph has this set to %d. Incompatible binary.\n",
			     file, num_b2, NUMBER_OF_BITFIELDS_IN_BINARY_KMER);
		      exit(1);
		    }
		  else if (num_c2 !=1)
		    {
		      printf("Input error. --colour_list requires a list of colour files, each one containing a list of single-colour binaries\n");
		      printf("This binary %s is not a single colour binary - it has %d colours.\n", file, num_c2);
		      exit(1);
		    }
		  count++;

		}
	    }
	  fclose(fp);




	}


      fclose(fptr);
    }


  //cleanup
  for (i=0; i<num_cols_in_list; i++)
    {
      free(list_cols[i]);
    }
  free(list_cols);

  return true;
}
