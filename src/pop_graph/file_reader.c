#include <binary_kmer.h>
#include <stdlib.h>
#include <stdio.h>
#include <dB_graph.h>
#include <seq.h>
#include <file_reader.h>
#include <global.h>
#include <string.h>

//takes a filename 
// this fle contains a list of filenames, each of these represents an individual (and contains a list of fasta for that individual).
int load_population_as_fasta(char* filename, dBGraph* db_graph)
{
  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
    printf("cannot open file:%s\n",filename);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
  }


  int MAX_FILENAME_LENGTH=300;
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
      //      if (people_so_far>NUMBER_OF_INDIVIDUALS_PER_POPULATION)
      //{
      //  printf("This filelist contains too many people for a single population, %d", people_so_far);
      //  exit(1);
      //}

      total_seq_loaded = total_seq_loaded + load_all_fasta_for_given_person_given_filename_of_file_listing_their_fasta_files(line, db_graph, people_so_far-1);

      //printf("Just loaded person number %d, and now have cumulative total of  %d bases\n", people_so_far-1, total_seq_loaded);
    }


  return total_seq_loaded;
  
}


//index tells you which person within a population it is
int load_all_fasta_for_given_person_given_filename_of_file_listing_their_fasta_files(char* f_name, dBGraph* db_graph, int index)
{
  FILE* fptr = fopen(f_name, "r");
  if (fptr == NULL)
    {
      //printf("cannot open person-specific fasta file:%s\n",f_name);
    exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
    }

  //file contains a list of fasta file names
  int MAX_FILENAME_LENGTH=300;
  char line[MAX_FILENAME_LENGTH+1];
  
  int total_seq_loaded=0;
  
  while(fgets(line,MAX_FILENAME_LENGTH, fptr) !=NULL)
    {

      //remove newline from endof line- replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';
      
      //open this fasta file
      FILE* temp_fptr = fopen(line, "r");
      if (temp_fptr == NULL)
	{
	  printf("cannot open  fasta file:%s\n",line);
	  exit(1); //TODO - prfer to print warning and skip file and reutnr an error code?
	}

      total_seq_loaded = total_seq_loaded + load_fasta_data_into_graph_for_specific_person_or_population(temp_fptr, db_graph, individual_edge_array, index);
    }

  return total_seq_loaded;

}



/*
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
	    }
	  
	  }
	  binary_kmer_free_kmers(&kmers);
	}
    }
  
  //fprintf(stderr, "Found this many bad reads:%d\n", count_bad_reads);

  return seq_length;
}

*/



//returns length of sequence loaded
int load_fasta_data_into_graph_for_specific_person_or_population(FILE* fp, dBGraph * db_graph, EdgeArrayType type, int index)
{

  //printf("ZAM start of  load_fasta_data_into_graph_for_specific_person_or_population\n");
  if (fp ==NULL)
    {
      printf("Do not give NUL pointer to load_fasta_data_into_graph_for_specific_person_or_population\n");
      exit(1);
    }

  Sequence* seq;
  int seq_length=0;
  int count_bad_reads=0;

  while ((seq = read_sequence_from_fasta(fp)))
    {
      if (DEBUG)
      	{
      printf ("\nZIM ZAM sequence %s\n",seq->seq);
	}

      KmerArray *kmers;
      int i;
      seq_length += seq->length;
      //printf("just got another line and so far in this file Seq length is %d\n", seq_length);
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
		printf("ZIMZAM kmer %i:  %s\n",i,binary_kmer_to_seq(kmers->bin_kmers[i],db_graph->kmer_size));
	      }
	  
	    if (i>0){
	      //never assume that previous pointer stays as we do reallocation !!!!!!
	      previous_node = hash_table_find(element_get_key(kmers->bin_kmers[i-1],db_graph->kmer_size),db_graph);
	      
	      if (previous_node == NULL){
		puts("file_reader: problem - kmer not found\n");
		exit(1);
	      }
	      previous_orientation = db_node_get_orientation(kmers->bin_kmers[i-1],previous_node, db_graph->kmer_size); 	      
	      db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size, type, index);	  	      
	      //printf("Add edge between %s and %s", binary_kmer_to_seq(previous_node->kmer, db_graph->kmer_size), binary_kmer_to_seq(current_node->kmer, db_graph->kmer_size) );
	    }
	  
	  }
	  binary_kmer_free_kmers(&kmers);
	}
    }
  
  //fprintf(stderr, "Found this many bad reads:%d\n", count_bad_reads);
  //printf("Seq length is %d\n", seq_length);
  return seq_length;
}


/*

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

*/
