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


  int count_file   = 0;
  long long total_length = 0; //total sequence length

  //Go through all the files, loading data into the graph
  while (!feof(fp_fnames)){

    //int count_bad_reads = 0;
    fscanf(fp_fnames, "%s\n", filename);
    
    int seq_length = 0;
    count_file++;

    total_length += load_population_as_fasta(filename, db_graph);

    fprintf(stderr,"\n%i file name:%s seq:%i total seq:%qd\n\n",count_file,filename,seq_length, total_length);
    
  }

  printf("print supernodes for each person in population 1\n");
  printf("print supernodes for person 1 population 1\n");

  //remember last two arguments are only used for unit tests, so here they are NULL, and 0.
  db_graph_traverse_specific_person_or_pop(&print_supernode_for_specific_person_or_pop,db_graph,individual_edge_array,0, NULL,0);
  db_graph_set_all_visited_nodes_to_status_none(db_graph);

  printf("print supernodes for person 2 population 1\n");
  db_graph_traverse_specific_person_or_pop(&print_supernode_for_specific_person_or_pop,db_graph,individual_edge_array,1, NULL,0);


  return 1;
}
