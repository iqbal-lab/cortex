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

int load_seq_data_into_graph_of_specific_person_or_pop(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length), long long * bad_reads, char quality_cut_off, int max_read_length, dBGraph * db_graph, EdgeArrayType type, int index);


int load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * bad_reads, int max_read_length, dBGraph* db_graph, EdgeArrayType type, int index)
{

  // printf("start to load %s\n", filename);
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  return load_seq_data_into_graph_of_specific_person_or_pop(fp,&read_sequence_from_fasta,bad_reads,  0 , max_read_length, db_graph, type, index);
}

int load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * bad_reads,  char quality_cut_off, int max_read_length, dBGraph* db_graph,
EdgeArrayType type, int index)
{
  //  printf("start to load %s\n", filename);
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  int seq_loaded = load_seq_data_into_graph_of_specific_person_or_pop(fp,&read_sequence_from_fastq, bad_reads, quality_cut_off, max_read_length, db_graph, type, index);

  //print mem status
  printf("Finished loading: %s, total bases loaded were: %d, bad reads: %lld, quality cutoff: %c, \nCurrent memory status:\n", filename, seq_loaded, *bad_reads, quality_cut_off);
  FILE* fmem=fopen("/proc/self/status", "r");
  char memline[500];
  while (fgets(memline,500,fmem) !=NULL){
    if (memline[0] == 'V' && memline[1] == 'm'){
      fprintf(stdout,"%s",memline);
    }
  }
  fclose(fmem);
  printf("************\n");

  return seq_loaded;
}



//returns length of sequence loaded
int load_seq_data_into_graph_of_specific_person_or_pop(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length), long long * bad_reads, char quality_cut_off, 
						       int max_read_length, dBGraph * db_graph, EdgeArrayType type, int index){

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
	  //printf("insertion. j is %d\n", j);
	  boolean found = false;
	  current_node = hash_table_find_or_insert(element_get_key(current_window->kmer[j],db_graph->kmer_size),&found,db_graph);	  	 
	  //if (!found){  //commented out - not counting new kmers any more
	  //  
	  //}

	  if (current_node==NULL)
	    {
	      printf("failed to find or insert node");
	      exit(1);
	    }
	  
	  //increment coverage, if not seen before in this read
	  if (!db_node_check_status(current_node, visited))
	    {
	      db_node_increment_coverage(current_node, type, index);
	      db_node_set_status(current_node, visited);
	    }

	  current_orientation = db_node_get_orientation(current_window->kmer[j],current_node, db_graph->kmer_size);
	  
	  if (DEBUG){
	    char kmer_seq[db_graph->kmer_size];
	    printf("kmer %i:  %s\n",i,binary_kmer_to_seq(current_window->kmer[j],db_graph->kmer_size,kmer_seq));
	  }
	  
	  if (j>0)
	    {
	      
	      if (previous_node == NULL)
		{
		  puts("file_reader: problem - kmer not found\n");
		  exit(1);
		}
	      else
		{
		  //previous_orientation = db_node_get_orientation(current_window->kmer[j-1],previous_node, db_graph->kmer_size); 	      
		  db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size, type, index);	  	      
		}
	    }
	  
	  previous_node = current_node;
          previous_orientation = current_orientation;

	  
	}
      }

      //here - go through all the windows again, find each node, and set it to unvisited.
      for(i=0;i<windows->nwindows;i++)
	{ 
	  //for each window
	  KmerSlidingWindow * current_window = &(windows->window[i]);
	  
	  for(j=0;j<current_window->nkmers;j++)
	    { //for each kmer in window
	      current_node = hash_table_find(element_get_key(current_window->kmer[j],db_graph->kmer_size),db_graph);	  	 
	      if (current_node==NULL)
		{
		  printf("i is %d and j is %d Impossible- must be able to find all kmers in these windows as already inserted", i, j);
		  exit(1);
		}
	      db_node_set_status(current_node, none);
	    }
	}
    }
  }
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
    
  return seq_length;    
}



//returns length of sequence loaded
int load_ref_overlap_data_into_graph_of_specific_person_or_pop(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length), long long * bad_reads, char quality_cut_off, 
						       int max_read_length, dBGraph * db_graph, int which_chromosome){

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

      
      Orientation current_orientation;
      
      for(i=0;i<windows->nwindows;i++){ //for each window
	KmerSlidingWindow * current_window = &(windows->window[i]);
	
	for(j=0;j<current_window->nkmers;j++){ //for each kmer in window

	  current_node = hash_table_find(element_get_key(current_window->kmer[j],db_graph->kmer_size),db_graph);	  	 

	  if (current_node){
	    current_orientation = db_node_get_orientation(current_window->kmer[j],current_node, db_graph->kmer_size);	    
	    db_node_mark_chromosome_overlap(current_node, which_chromosome, current_orientation);
	  }
	  
	}
	
      }
    }
    
  }
  
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);


  //print mem status
  //FILE* fmem=fopen("/proc/self/status", "r");
  //char memline[500];
  //while (fgets(memline,500,fmem) !=NULL){
  //  if (memline[0] == 'V' && memline[1] == 'm'){
  //    fprintf(stderr,"%s",memline);
  //  }
  // }
  //fclose(fmem);
  //fprintf(stderr,"************\n");


  return seq_length;    
}




//takes a filename 
// this file contains a list of filenames, each of these represents an individual (and contains a list of fasta for that individual).
long long load_population_as_fasta(char* filename, long long* bad_reads, dBGraph* db_graph)
{

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
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
  
  while(fgets(line,MAX_FILENAME_LENGTH, fptr) !=NULL)
    {

      //remove newline from endof line- replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';
      
      total_seq_loaded = total_seq_loaded + 
	load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(line, bad_reads, MAX_READ_LENGTH, db_graph, individual_edge_array, index);

    }

  return total_seq_loaded;

}




//takes a filename 
// this file contains a list of filenames, each of these represents an individual (and contains a list of fasta for that individual).
long long load_population_as_fastq(char* filename, long long* bad_reads, char quality_cutoff, dBGraph* db_graph)
{

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
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
  
  int total_seq_loaded=0;
  
  while(fgets(line,MAX_FILENAME_LENGTH, fptr) !=NULL)
    {

      //remove newline from endof line- replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';
      
      total_seq_loaded = total_seq_loaded + 
	load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(line, bad_reads, quality_cutoff,  MAX_READ_LENGTH, db_graph, individual_edge_array, index);

    }

  return total_seq_loaded;

}


//assumes you already have loaded your reads from your individuals, so have a full graph.
//this just marks each node with how it overlaps with this chromosome
int load_chromosome_overlap_data(char* f_name,  dBGraph* db_graph, int which_chromosome)
{
  
  int max_read_length = 500;

  FILE* fptr = fopen(f_name, "r");
  if (fptr == NULL)
    {
      //printf("cannot open chromosome fasta file:%s\n",f_name);
      exit(1);
    }
  //else
  //  {
  //    printf("opened file %s to load overlap data\n", f_name);
  //  }

  long long bad_reads;
  int total_seq_loaded=0;

  total_seq_loaded = total_seq_loaded +
    load_ref_overlap_data_into_graph_of_specific_person_or_pop(fptr, &read_sequence_from_fasta, &bad_reads,  0 , max_read_length, db_graph, which_chromosome);


  return total_seq_loaded;

}


//loading a binary as dumped by sv_trio. So contains information on the whole trio.
//returns number of kmers loaded
int load_sv_trio_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph)
{
  FILE* fp_bin = fopen(filename, "r");
  int seq_length = 0;
  dBNode node_from_file;
  boolean found;
  int count=0;

  if (fp_bin == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }
  
  //Go through all the entries in the binary file
  while (db_node_read_sv_trio_binary(fp_bin,db_graph->kmer_size,&node_from_file)){
    count++;
    
    //if (count % 100000000 == 0 ){
    // printf("loaded %i\n",count);

    //}
   
    dBNode * current_node  = NULL;
    current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size),&found,db_graph);
    
    seq_length+=db_graph->kmer_size;
   
    int i;
    for (i=0; i<NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
      {
	add_edges(current_node,individual_edge_array, i, get_edge_copy(node_from_file, individual_edge_array, i));
	db_node_update_coverage(current_node, individual_edge_array, i, db_node_get_coverage_as_unsigned_char(&node_from_file, individual_edge_array,i));
      }
    for (i=0; i<6; i++)
      {
	current_node->chrom_xs[i] |= node_from_file.chrom_xs[i];
      }
  }
  
  fclose(fp_bin);
  return seq_length;
}




//reads a binary as dumped by graph (not sv_trio)
int load_individual_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph, EdgeArrayType type, int index)
{

  FILE* fp_bin = fopen(filename, "r");
  int seq_length = 0;
  dBNode node_from_file;
  boolean found;
  int count=0;

  if (fp_bin == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }
  
  //Go through all the entries in the binary file
  while (db_node_read_graph_binary(fp_bin,db_graph->kmer_size,&node_from_file, type, index))
    {
      count++;

      //if (count % 100000000 == 0 ){
      // printf("loaded %i\n",count);
      //}
      
      dBNode * current_node  = NULL;
      current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size),&found,db_graph);
      
      seq_length+=db_graph->kmer_size;//todo - maybe only increment if had to insert, not if was already in graph?
      
      add_edges(current_node,individual_edge_array, index, get_edge_copy(node_from_file, individual_edge_array, index));
      db_node_update_coverage(current_node, individual_edge_array, index, db_node_get_coverage_as_unsigned_char(&node_from_file, individual_edge_array,index) );
      int i;
      for (i=0; i<6; i++)
	{
	  current_node->chrom_xs[i] |= node_from_file.chrom_xs[i];
	}
    }
  
  fclose(fp_bin);
  return seq_length;

}
