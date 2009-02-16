
#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>

int main(int argc, char **argv){

  FILE *fp_bin;

  int hash_key_bits;
  dBGraph * db_graph = NULL;
  dBNode node_from_file;
  short kmer_size;
  int clip_limit;
  long long kmers_clipped = 0;
  //long long count_kmers   = 0;
  boolean found;

  void tip_clipping(dBNode * db_node){
    kmers_clipped += db_graph_clip_tip(db_node,clip_limit,db_graph);
  }


  //command line arguments 
  fp_bin = fopen(argv[1], "r");    //open bin file
  kmer_size        = atoi(argv[2]); 
  hash_key_bits    = atoi(argv[3]); //number of buckets: 2^hash_key_bits
  clip_limit       = atoi(argv[4]);
  DEBUG            = atoi(argv[7]);
  
  fprintf(stderr,"Kmer size: %d hash_table_size (%d bits): %d\n",kmer_size,hash_key_bits,1 << hash_key_bits);
 
  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(1 << hash_key_bits,kmer_size);
  fprintf(stderr,"table created: %d\n",1 << hash_key_bits);

  
  //Go through all the entries in the binary file
  while (db_node_read_binary(fp_bin,kmer_size,&node_from_file)){
    dBNode * current_node  = NULL;

    current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file),kmer_size),&found,db_graph);
    
    db_node_set_edges(current_node,db_node_get_edges(&node_from_file));
  }
					     //    fprintf(stderr,"\n%i kmers: %qd file name:%s bad reads: %qd seq:%i total seq:%qd\n\n",count_file,count_kmers,filename,bad_reads,seq_length, total_length);

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
  
  if (clip_limit>0){
    printf("clip tips\n");
    hash_table_traverse(&tip_clipping,db_graph);
    printf("\nkmers clipped %qd\n\n",kmers_clipped);
  }
  
  fclose(fp_bin);
  return 1;
}
