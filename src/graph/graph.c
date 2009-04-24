#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>

int main(int argc, char **argv){

  FILE *fp_fnames;
  char filename[1000];
  int hash_key_bits;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int bucket_size;
  int action; //0 dump graph - 1 call SNPs
  char* binfilename;

  FILE * fout; //binary output

   //command line arguments 
  fp_fnames= fopen(argv[1], "r");    //open file of file names
  kmer_size        = atoi(argv[2]); 
  hash_key_bits    = atoi(argv[3]);  //number of buckets: 2^hash_key_bits
  bucket_size      = atoi(argv[4]);
  action           = atoi(argv[5]);
  DEBUG            = atoi(argv[6]);
  binfilename      = argv[7];


  fprintf(stderr,"Input file of filenames: %s - action: %i\n",argv[1],action);
  fprintf(stderr,"Kmer size: %d hash_table_size (%d bits): %d - bucket size: %d - total size: %qd\n",kmer_size,hash_key_bits,1 << hash_key_bits, bucket_size, ((long long) 1<<hash_key_bits)*bucket_size);

  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(hash_key_bits,bucket_size, 10,kmer_size);
  fprintf(stderr,"table created: %d\n",1 << hash_key_bits);

  int count_file   = 0;
  long long total_length = 0; //total sequence length

  //Go through all the files, loading data into the graph
  while (!feof(fp_fnames)){

    fscanf(fp_fnames, "%s\n", filename);
    
    long long seq_length = 0;
    count_file++;

    seq_length += load_binary_data_from_filename_into_graph(filename,db_graph);

    total_length += seq_length;
    
    fprintf(stderr,"\n%i kmers: %qd file name:%s seq:%qd total seq:%qd\n\n",count_file,hash_table_get_unique_kmers(db_graph),filename,seq_length, total_length);

    //print mem status
    FILE* fmem=fopen("/proc/self/status", "r");
    char line[500];
    while (fgets(line,500,fmem) !=NULL){
      if (line[0] == 'V' && line[1] == 'm'){
	fprintf(stderr,"%s",line);
      }
    }
    fclose(fmem);
    fprintf(stderr,"************\n");
  }  
    
  
 

  
  switch (action){
  case 0 :
    printf("dumping graph %s\n",argv[7]);
    fout= fopen(argv[7], "w"); 

    //routine to dump graph
    void print_node_binary(dBNode * node){
    db_node_print_binary(fout,node);
    }

    hash_table_traverse(&print_node_binary,db_graph); 
    break;

  case 1 :
    printf("call SNPs\n");
    db_graph_detect_snps(db_graph);
    break;

  case 2:
    printf("clip tips\n");
    db_graph_clip_tips(db_graph);
    break;

  case 3:
    printf("print supernodes\n");
    db_graph_print_supernodes(db_graph); 
    break;

  case 4:
    printf("clip tips\n");
    db_graph_clip_tips(db_graph);
    printf("print supernodes\n");
    db_graph_print_supernodes(db_graph); 
    break;

  case 5:
    printf("detect SNPs\n");
    db_graph_detect_snps(db_graph);
    break;
    
  case 6:
    printf("count kmers\n");
    db_graph_print_coverage(db_graph);
    break;
    
  case 7:
    printf("Print supernodes containing entirely novel(non-reference) sequence\n");
    read_all_ref_chromosomes_and_mark_graph(db_graph);
    int min_covg=2; //demand coverage of at least 2

    //todo - add condition to make sure ignore pruned

    db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(db_graph, &db_node_check_status_is_not_exists_in_reference,min_covg, false, NULL, 0);

    hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph); //cleanup - the printing process has set all printed nodes to visited
                                                                                                            //        - or to visited_and_exists_in_reference
  }

  return 0;
}
