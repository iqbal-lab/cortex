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

int load_seq_data_into_graph(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length), long long * count_kmers, long long * bad_reads, char qualiy_cut_off, int max_read_length, dBGraph * db_graph);


int load_fasta_data_from_filename_into_graph(char* filename, long long * count_kmers, long long * bad_reads, int max_read_length, dBGraph* db_graph)
{
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  int ret = load_seq_data_into_graph(fp,&read_sequence_from_fasta, count_kmers,bad_reads,  0 , max_read_length, db_graph);

  fclose(fp);

  return ret;
}

int load_fastq_data_from_filename_into_graph(char* filename, long long * count_kmers, long long * bad_reads,  char quality_cut_off, int max_read_length, dBGraph* db_graph)
{
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  int ret =  load_seq_data_into_graph(fp,&read_sequence_from_fastq, count_kmers, bad_reads, quality_cut_off, max_read_length, db_graph);
  fclose(fp);

  return ret;
}




//returns length of sequence loaded
int load_seq_data_into_graph(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length), long long * count_kmers, long long * bad_reads, char quality_cut_off, int max_read_length, dBGraph * db_graph){

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
  
  int entry_length;

  while ((entry_length = file_reader(fp,seq,max_read_length))){

    if (DEBUG){
      printf ("\nsequence %s\n",seq->seq);
    }
    
    int i,j;
    seq_length += entry_length;
    
    int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,quality_cut_off,db_graph->kmer_size,windows,max_windows, max_kmers);
    
    if (nkmers == 0) {
      (*bad_reads)++;
    }
    else {
      Element * current_node  = NULL;
      Element * previous_node = NULL;
      
      Orientation current_orientation,previous_orientation;
      
      for(i=0;i<windows->nwindows;i++){ //for each window
	KmerSlidingWindow * current_window = &(windows->window[i]);
	
	for(j=0;j<current_window->nkmers;j++){ //for each kmer in window
	  boolean found = false;
	  current_node = hash_table_find_or_insert(element_get_key(current_window->kmer[j],db_graph->kmer_size),&found,db_graph);	  	 
	  if (!found){
	    (*count_kmers)++;
	  }
	  
	  current_orientation = db_node_get_orientation(current_window->kmer[j],current_node, db_graph->kmer_size);
	  
	  if (DEBUG){
	    char kmer_seq[db_graph->kmer_size];
	    printf("kmer %i:  %s\n",i,binary_kmer_to_seq(current_window->kmer[j],db_graph->kmer_size,kmer_seq));
	  }
	  
	  if (j>0){
	    //never assume that previous pointer stays as we do reallocation !!!!!!
	    previous_node = hash_table_find(element_get_key(current_window->kmer[j-1],db_graph->kmer_size),db_graph);
	    
	    if (previous_node == NULL){
	      puts("file_reader: problem - kmer not found\n");
	      exit(1);
	    }
	    previous_orientation = db_node_get_orientation(current_window->kmer[j-1],previous_node, db_graph->kmer_size); 	      
	    db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size);	  	      
	  }
	  
	}
      }
     
    }
  }
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
  return seq_length;    
}


//returns number of sequence loaded (ie all the kmers concatenated)
//count_kmers returns the number of new kmers
int load_binary_data_from_filename_into_graph(char* filename, long long * count_kmers, dBGraph* db_graph)
{
  FILE* fp_bin = fopen(filename, "r");
  int seq_length = 0;
  dBNode node_from_file;
  boolean found;

  if (fp_bin == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  //Go through all the entries in the binary file
  while (db_node_read_binary(fp_bin,db_graph->kmer_size,&node_from_file)){
    dBNode * current_node  = NULL;
    current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size),&found,db_graph);
    if (!found){
      (*count_kmers)++;
    }
    
    seq_length+=db_graph->kmer_size;
   
    db_node_set_edges(current_node,db_node_get_edges(&node_from_file));
  }
  
  fclose(fp_bin);
  return seq_length;
}

