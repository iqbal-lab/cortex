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

int load_seq_data_into_graph_of_specific_person_or_pop(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
						       long long * bad_reads, char quality_cut_off, int max_read_length, dBGraph * db_graph, EdgeArrayType type, int index);


//this routine supports big fasta entries (chromosome length for example)
int load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * bad_reads, int max_chunk_length, dBGraph* db_graph, EdgeArrayType type, int index)
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
    fprintf(stderr," load_fasta_data_from_filename_into_graph_of_specific_person_or_pop cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  long long ret =  load_seq_data_into_graph_of_specific_person_or_pop(fp,&file_reader,bad_reads,0,max_chunk_length,db_graph, type, index);
  fclose(fp);


  return ret;

}

int load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * bad_reads,  char quality_cut_off, int max_read_length, dBGraph* db_graph,
								       EdgeArrayType type, int index)
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
    fprintf(stderr,"load_fastq_data_from_filename_into_graph_of_specific_person_or_pop cannot open file:%s\n",filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }

  long long ret =  load_seq_data_into_graph_of_specific_person_or_pop(fp,&file_reader,bad_reads,quality_cut_off,max_read_length,db_graph, type, index);
  fclose(fp);

  return ret;


}



//returns length of sequence loaded
int load_seq_data_into_graph_of_specific_person_or_pop(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
						       long long * bad_reads, char quality_cut_off, 
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


  while ((entry_length = file_reader(fp,seq,max_read_length, full_entry, &full_entry))){

    if (DEBUG){
      printf("entry %s %i %i %s\n",seq->name,seq->start,seq->end,prev_full_entry?"not connecting to previous read in same fasta entry":"connecting to previous read in same fasta entry");
    }

    int i,j;
    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));

    //printf("Length %qd %i %s %i %i\n",seq_length,entry_length, seq->name, seq->start, seq->end);

    int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,quality_cut_off,db_graph->kmer_size,windows,max_windows, max_kmers);

    //printf("number kmers:%i %i\n",nkmers,windows->nwindows);
    
    if (nkmers == 0) 
      {
	(*bad_reads)++;
      }
    else 
      {
	
	Element * current_node  = NULL;
	Element * previous_node  = NULL;
	Orientation current_orientation;
	Orientation previous_orientation;
	
	for(i=0;i<windows->nwindows;i++)
	  { //for each window
	    KmerSlidingWindow * current_window = &(windows->window[i]);
	    
	    for(j=0;j<current_window->nkmers;j++) //for each kmer in window
	      {
		boolean found = false;
		current_node = hash_table_find_or_insert(element_get_key(current_window->kmer[j],db_graph->kmer_size),&found,db_graph);	  
		if (current_node == NULL)
		  {
		    fputs("file_reader: problem - current kmer not found\n",stderr);
		    exit(1);
		  }
		
		if (! (i==0 && j==0 && prev_full_entry == false && current_node == previous_node)) //otherwise is the same old last entry
		  {
		    db_node_increment_coverage(current_node, type, index);
		  }

		current_orientation = db_node_get_orientation(current_window->kmer[j],current_node, db_graph->kmer_size);
		
		if (DEBUG)
		  {
		    char kmer_seq[db_graph->kmer_size];
		    char kmer_seq2[db_graph->kmer_size];
		    printf("kmer i:%i j:%i:  %s %s %i\n",i,j,binary_kmer_to_seq(current_window->kmer[j],db_graph->kmer_size,kmer_seq),
			   binary_kmer_to_seq(binary_kmer_reverse_complement(current_window->kmer[j],db_graph->kmer_size),db_graph->kmer_size,kmer_seq2),
			   db_node_get_coverage(current_node, type, index));
		  }
		
		if (j>0)
		  {
		    
		    if (previous_node == NULL)
		      {
			printf("PROBLEM i:%i j:%i seq_length:%qd nkmers:%i prev_full_entry:%s\n",i,j,seq_length,nkmers,prev_full_entry == true ? "true" : "false");
			fputs("file_reader: problem - prev kmer not found\n",stderr);
			exit(1);
		      }
		    else
		      {
			db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size, type, index);	  	  
		      }
		  }
		previous_node = current_node;
		previous_orientation = current_orientation;
		
	      }
	  }
      }//matches else..
    
    if (full_entry == false)
      {
	shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
      }
    
    prev_full_entry = full_entry;
    

  }

  

  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
  
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

  if (DEBUG){
    printf ("\nsequence returned from read_sequence_from_fasta to load_seq_into_array is %s - kmer size: %i - number of bases loaded inc the preassigned ones at start f seq  length: %i \n",seq->seq,db_graph->kmer_size, chunk_length);
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

  //if num_kmers=0 (ie have hit end of file, then kmer_window->nkmers=0), so will skip this
  for(j=0;j<kmer_window->nkmers;j++)
    { //for each kmer in window

      if (kmer_window->kmer[j]==~0) //encoding as  1's in all 64 bits, a non-kmer - ie anything that would have had an N in it
	{
	  //corresponds to a kmer that contains an N
	  path_nodes[offset+j]        =NULL;
	  path_orientations[offset+j] =forward;
	  path_labels[offset+j]       =Undefined;
	  //path_string[offset+j]       ='N';
	  previous_node=NULL;
	  continue;
	}
      
      boolean found = false;
      current_node = hash_table_find_or_insert(element_get_key(kmer_window->kmer[j],db_graph->kmer_size),&found,db_graph);	  

      if (current_node == NULL){
	fputs("Problem in load_seq_into_array - current kmer not found\n",stderr);
	exit(1);
      }
      
      current_orientation = db_node_get_orientation(kmer_window->kmer[j],current_node, db_graph->kmer_size);

      path_nodes[offset+j]        = current_node;
      path_orientations[offset+j] = current_orientation;


      if (DEBUG){
	char kmer_seq[db_graph->kmer_size];
	printf("j=%d, Current node kmer is  %s\n",j, binary_kmer_to_seq(kmer_window->kmer[j],db_graph->kmer_size,kmer_seq));
	if (current_orientation==forward)
	  printf("Current orientation is forward\n");
	else
	  printf("Current orientation is reverse");
      }
      
      if (j>0)
	{
	  if (previous_node == NULL)
	    {
	      path_labels[offset+j]=Undefined;
	      //path_string[offset+j]=='N';
	      //path_string[offset+j]='N';
	      //path_string[offset+j+1]='\0';
	    }
	  else
	    {
	      BinaryKmer previous_k, current_k; 
	      char seq1[kmer_size];
	      char seq2[kmer_size];
	      
	      previous_k = previous_node->kmer;
	      current_k  = current_node->kmer;
	      
	      if (previous_orientation == reverse){
		previous_k = binary_kmer_reverse_complement(previous_k,kmer_size);
	      }
	      
	      if (current_orientation == reverse){
		current_k = binary_kmer_reverse_complement(current_k,kmer_size);
	      }
    
	      
	      if (DEBUG){
		printf("Found edge %s -%c-> %s\n",binary_kmer_to_seq(previous_k,kmer_size,seq1),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(current_k)),binary_kmer_to_seq(current_k,kmer_size,seq2));
	      }
	      path_labels[offset+j]=binary_kmer_get_last_nucleotide(current_k);
	      //path_string[offset+j]=binary_nucleotide_to_char(path_labels[offset+j]);
	      //path_string[offset+j+1]='\0';
	      
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


/*
//returns length of sequence loaded
int load_ref_overlap_data_into_graph_of_specific_person_or_pop(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry, int offset), 
							       long long * bad_reads, char quality_cut_off,  int max_read_length, dBGraph * db_graph, int which_chromosome)
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
  
  int entry_length;
  boolean full_entry = true;
  boolean prev_full_entry = true;


  while ((entry_length = file_reader(fp,seq,max_read_length, full_entry, &full_entry))){

    if (DEBUG){
      printf ("\nsequence %s\n",seq->seq);
    }
    
    int i,j;
    seq_length  += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));

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

    if (full_entry == false)
      {
        shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
      }
    
    prev_full_entry = full_entry;
    
    
  }
  
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);


  return seq_length;    
}

*/


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
  
  while(fgets(line,MAX_FILENAME_LENGTH, fptr) !=NULL)
    {

      //remove newline from endof line- replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';
      
      total_seq_loaded = total_seq_loaded + 
	load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(line, bad_reads, MAX_READ_LENGTH, db_graph, individual_edge_array, index);

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

  fclose(fptr);
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
    printf("load_sv_trio_binary_data_from_filename_into_graph cannot open file:%s\n",filename);
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
	db_node_update_coverage(current_node, individual_edge_array, i, db_node_get_coverage_as_short(&node_from_file, individual_edge_array,i));
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
    printf("load_individual_binary_data_from_filename_into_graph cannot open file:%s\n",filename);
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
      db_node_update_coverage(current_node, individual_edge_array, index, db_node_get_coverage_as_short(&node_from_file, individual_edge_array,index) );
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
	load_individual_binary_data_from_filename_into_graph(line, db_graph, individual_edge_array, index);

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
