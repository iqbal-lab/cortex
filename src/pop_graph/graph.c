#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <dB_graph_population.h>

int main(int argc, char **argv){


  FILE *fp_fnames;
  char filename[100];
  int hash_key_bits;
  dBGraph * db_graph = NULL; 
  short kmer_size;


  //command line arguments 
  fp_fnames= fopen(argv[1], "r");    //open file of which gives one filename per population. Each of these files gives one filename per individual, and each of those gives list of files.
  kmer_size        = atoi(argv[2]);  //global variable defined in element.h
  hash_key_bits    = atoi(argv[3]); //number of buckets: 2^hash_key_bits
  DEBUG            = atoi(argv[4]);

  fprintf(stderr,"Kmer size: %d hash_table_size (%d bits): %d\n",kmer_size,hash_key_bits,1 << hash_key_bits);

  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(1 << hash_key_bits,kmer_size);
  fprintf(stderr,"table created: %d\n",1 << hash_key_bits);

  long count_files = 0;
  long long bad_reads =0;
  long long total_kmers =0;
  long long total_length = 0; //total sequence length

  //Go through all the files, loading data into the graph
  while (!feof(fp_fnames)){

    fscanf(fp_fnames, "%s\n", filename);
    count_files++;
    total_length += load_population_as_fasta(filename, &total_kmers, &bad_reads, db_graph);
    fprintf(stderr,"\n%ld file name:%s total seq:%lld bad reads: %lld total kmers:%lld\n\n",count_files,filename,total_length, bad_reads, total_kmers);
    
  }


  //remember last two arguments are only used for unit tests, so here they are NULL, and 0.
  long file_count0=0;
  db_graph_traverse_specific_person_or_pop(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop,db_graph,file_count0, individual_edge_array,0, false,NULL,0);
  db_graph_set_all_visited_nodes_to_status_none(db_graph);

  printf("print supernodes for person 2 population 1\n");
  long file_count1=0;
  db_graph_traverse_specific_person_or_pop(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop,db_graph,file_count1, individual_edge_array,1, false, NULL,0);


  int NUM_PEOPLE=2;
  int stats[NUM_PEOPLE+1];
  int* array[NUM_PEOPLE+1];

  int j;
  for (j=0; j<=NUM_PEOPLE; j++)
    {
      stats[j]=0;
      array[j]=&stats[j];
    }

  hash_table_traverse_2(&find_out_how_many_individuals_share_this_node_and_add_to_statistics , db_graph, array, NUM_PEOPLE);

  printf("Results are:\n");
  for (j=0; j<=NUM_PEOPLE; j++)
    {
      printf("Number of kmers shared by %d people is %d\n", j, stats[j]);
    }

  

  return 1;
}
