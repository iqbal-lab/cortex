#include <binary_kmer.h>
#include <stdlib.h>
#include <stdio.h>
#include <dB_graph.h>
#include <seq.h>
#include <file_reader.h>

int load_fasta_data_from_filename_into_graph_efficient(char* filename, dBGraph* db_graph, int* number_of_reallocs)
{
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
  }

  return load_fasta_data_into_graph_efficient(fp, db_graph, number_of_reallocs);
}

int load_fasta_data_from_filename_into_graph(char* filename, dBGraph* db_graph)
{
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
  }

  return load_fasta_data_into_graph(fp, db_graph);
}


int load_fastq_data_from_filename_into_graph(char* filename, dBGraph* db_graph)
{
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
  }

  return load_fastq_data_into_graph(fp, db_graph);
}


//returns length of sequence loaded
int load_fasta_data_into_graph(FILE* fp, dBGraph * db_graph)
{
  Sequence* seq;
  int seq_length=0;
  int count_bad_reads=0;

  while ((seq = read_sequence_from_fasta(fp)))
    {
      if (DEBUG)
	{
	  printf ("\nsequence %s\n",seq->seq);
	}

      KmerArray *kmers;
      int i;
      seq_length += seq->length;

      kmers = get_binary_kmers_from_sequence(seq->seq,seq->length,db_graph->kmer_size);
      free_sequence(&seq);

      if (kmers == NULL)
	{
	  count_bad_reads++;
	}

      else
	{
	  Element * current_node  = NULL;
	  Element * previous_node = NULL;
 	
	  Orientation current_orientation,previous_orientation;
	
	  for(i=0;i<kmers->nkmers;i++){	   
     	    current_node = hash_table_find_or_insert(element_get_key(kmers->bin_kmers[i],db_graph->kmer_size),db_graph);	  	  
	    current_orientation = db_node_get_orientation(kmers->bin_kmers[i],current_node, db_graph->kmer_size);
	    
	    if (DEBUG)
	      {
		printf("kmer %i:  %s\n",i,binary_kmer_to_seq(kmers->bin_kmers[i],db_graph->kmer_size));
	      }
	  
	    if (i>0){
	      //never assume that previous pointer stays as we do reallocation !!!!!!
	      previous_node = hash_table_find(element_get_key(kmers->bin_kmers[i-1],db_graph->kmer_size),db_graph);
	      
	      if (previous_node == NULL){
		puts("file_reader: problem - kmer not found\n");
		exit(1);
	      }
	      previous_orientation = db_node_get_orientation(kmers->bin_kmers[i-1],previous_node, db_graph->kmer_size); 	      
	      db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size);	  	      
	      printf("Add edge from %s to %s\n", binary_kmer_to_seq(previous_node->kmer, db_graph->kmer_size), binary_kmer_to_seq(current_node->kmer, db_graph->kmer_size));
	    }
	  
	  }
	  binary_kmer_free_kmers(&kmers);
	}
    }
  
  //fprintf(stderr, "Found this many bad reads:%d\n", count_bad_reads);

  return seq_length;
}


//returns length of sequence loaded
int load_fasta_data_into_graph_efficient(FILE* fp, dBGraph * db_graph, int* number_of_reallocs)
{

  int longest_expected_read_length = 3000;

  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL)
    {
      printf("Out of memory allocating a Sequence object");
      exit(1);
    }
  //ZAM todo - replace woth calloc so zero initialised.
  seq->seq = malloc(sizeof(char)*longest_expected_read_length ) ;
  if (seq->seq == NULL)
    {
      printf("Out of memory allocating a seq->seq object");
      exit(1);
    }

  seq->max=longest_expected_read_length;
  seq->qual=NULL;
  seq->length=0;
  

  KmerArray * kmers = malloc(sizeof(KmerArray));  
  if (kmers == NULL){
    fputs("Out of memory trying to allocate a Kmer",stderr);
    exit(1);
  }


  kmers->bin_kmers = malloc(sizeof(BinaryKmer) * (longest_expected_read_length - db_graph->kmer_size + 1 ));
  if (kmers->bin_kmers == NULL){
    fputs("Out of memory trying to allocate a bin Kmers",stderr);
    exit(1);
  }



  int seq_length=0;
  int count_bad_reads=0;
  int total_reallocs=0;
  int did_a_realloc_occur_in_this_call=0;
  while ((read_sequence_from_fasta_efficient(fp,seq, &did_a_realloc_occur_in_this_call)))
    {
      total_reallocs=total_reallocs+did_a_realloc_occur_in_this_call;

      if (DEBUG)
	{
	  printf ("\nsequence %s\n",seq->seq);
	}

      //KmerArray *kmers;
      int i;
      seq_length += seq->length;

      int ret  = get_binary_kmers_from_sequence_efficient(seq->seq,seq->length,db_graph->kmer_size, kmers);
      seq->name=0;
      seq->length=0;
      //*(seq->seq)=0;
      //seq->max=0;
      //seq->qual=0;

      if (ret == 0)
	{
	  count_bad_reads++;
	}

      else
	{
	  Element * current_node  = NULL;
	  Element * previous_node = NULL;
 	
	  Orientation current_orientation,previous_orientation;
	
	  for(i=0;i<kmers->nkmers;i++){	   
     	    current_node = hash_table_find_or_insert(element_get_key(kmers->bin_kmers[i],db_graph->kmer_size),db_graph);	  	  
	    current_orientation = db_node_get_orientation(kmers->bin_kmers[i],current_node, db_graph->kmer_size);
	    
	    if (DEBUG)
	      {
		printf("kmer %i:  %s\n",i,binary_kmer_to_seq(kmers->bin_kmers[i],db_graph->kmer_size));
	      }
	  
	    if (i>0){
	      //never assume that previous pointer stays as we do reallocation !!!!!!
	      previous_node = hash_table_find(element_get_key(kmers->bin_kmers[i-1],db_graph->kmer_size),db_graph);
	      
	      if (previous_node == NULL){
		puts("file_reader: problem - kmer not found\n");
		exit(1);
	      }
	      previous_orientation = db_node_get_orientation(kmers->bin_kmers[i-1],previous_node, db_graph->kmer_size); 	      
	      db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size);	  	      
	    }
	  
	  }
	  //*(kmers->bin_kmers)=0;
	}
    }
  
  //fprintf(stderr, "Found this many bad reads:%d\n", count_bad_reads);

  free(seq->seq);
  free(seq);
  free(kmers->bin_kmers);
  free(kmers);

  *number_of_reallocs = total_reallocs;
  return seq_length;
}




//returns length of sequence loaded
int load_fastq_data_into_graph(FILE* fp, dBGraph * db_graph)
{
  Sequence* seq;
  int seq_length=0;
  int count_bad_reads=0;

  while ((seq = read_sequence_from_fastq(fp)))
    {
      if (DEBUG)
	{
	  printf ("\nsequence %s\n",seq->seq);
	}

      KmerArray *kmers;
      int i;
      seq_length += seq->length;

      kmers = get_binary_kmers_from_sequence(seq->seq,seq->length,db_graph->kmer_size);
      free_sequence(&seq);

      if (kmers == NULL)
	{
	  count_bad_reads++;
	}

      else
	{
	  Element * current_node  = NULL;
	  Element * previous_node = NULL;
 	
	  Orientation current_orientation,previous_orientation;
	
	  for(i=0;i<kmers->nkmers;i++){	   
     	    current_node = hash_table_find_or_insert(element_get_key(kmers->bin_kmers[i],db_graph->kmer_size),db_graph);	  	  
	    current_orientation = db_node_get_orientation(kmers->bin_kmers[i],current_node, db_graph->kmer_size);
	    
	    if (DEBUG)
	      {
		printf("kmer %i:  %s\n",i,binary_kmer_to_seq(kmers->bin_kmers[i],db_graph->kmer_size));
	      }
	  
	    if (i>0){
	      //never assume that previous pointer stays as we do reallocation !!!!!!
	      previous_node = hash_table_find(element_get_key(kmers->bin_kmers[i-1],db_graph->kmer_size),db_graph);
	      
	      if (previous_node == NULL){
		puts("file_reader: problem - kmer not found\n");
		exit(1);
	      }
	      previous_orientation = db_node_get_orientation(kmers->bin_kmers[i-1],previous_node, db_graph->kmer_size); 	      
	      db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size);	  	      
	    }
	  
	  }
	  binary_kmer_free_kmers(&kmers);
	}
    }
  
  fprintf(stderr, "Found this many bad reads:%d\n", count_bad_reads);

  return seq_length;
}


