
#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>

int main(int argc, char **argv){

  FILE *fp_fnames;
  char filename[100];
  int hash_key_bits;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int clip_limit;
  long long kmers_clipped = 0;
  long long count_kmers = 0;

  FILE * fout;
  long count=0;
  
  void print_supernode(dBNode * node){

    char filename [50];
    if (count % 100000000 == 0){
    int index = count / 100000000;

    if (count !=0){
      fclose(fout);
    }

    sprintf(filename,"out_nodes_cl_%i_%i_%i",kmer_size,clip_limit,index);
    fprintf(stderr,"opening file %s\n",filename);
    fout = fopen(filename,"w");
    }
    
    count++;
    db_graph_print_supernode(fout,node,db_graph);
  }


  void tip_clipping(dBNode * node){

    kmers_clipped += db_graph_clip_tip(node,clip_limit,db_graph);
    
  }


  //command line arguments 
  fp_fnames= fopen(argv[1], "r");    //open file of file names
  kmer_size        = atoi(argv[2]);  //global variable defined in element.h
 
  hash_key_bits    = atoi(argv[3]); //number of buckets: 2^hash_key_bits
  clip_limit       = atoi(argv[4]);
  DEBUG            = atoi(argv[5]);

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

    total_length += load_fasta_data_from_filename_into_graph(filename, &count_kmers,0, 5000, db_graph);

    fprintf(stderr,"\n%i kmers: %qd file name:%s seq:%i total seq:%qd\n\n",count_file,count_kmers,filename,seq_length, total_length);
    
  }

  printf("clip tips\n");
  hash_table_traverse(&tip_clipping,db_graph);

  printf("\nkmers clipped %qd\n\n",kmers_clipped);

  printf("print supernodes\n");
  hash_table_traverse(&print_supernode,db_graph);

  return 1;
}
