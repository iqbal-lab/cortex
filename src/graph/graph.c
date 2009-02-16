#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>

int main(int argc, char **argv){

  FILE *fp_fnames, *fp_file;
  char filename[300];

  int hash_key_bits;
  dBGraph * db_graph = NULL; 
  short kmer_size;

  FILE * fout;
  long count=0;

  void print_supernode(dBNode * node){

    char filename [50];
    if (count % 100000000 == 0){
    int index = count / 100000000;

    if (count !=0){
      fclose(fout);
    }

    sprintf(filename,"out_nodes_%i_%i",kmer_size,index);
    fprintf(stdout,"opening file %s\n",filename);
    fout = fopen(filename,"w");
    }
    
    count++;
    db_graph_print_supernode(fout,node,db_graph);
  }


  //command line arguments 
  fp_fnames= fopen(argv[1], "r");    //open file of file names
  kmer_size        = atoi(argv[2]);  //global variable defined in element.h
 
  hash_key_bits    = atoi(argv[3]); //number of buckets: 2^hash_key_bits
  DEBUG            = atoi(argv[4]);

  fprintf(stdout,"Kmer size: %d hash_table_size (%d bits): %d\n",kmer_size,hash_key_bits,1 << hash_key_bits);

  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(1 << hash_key_bits,kmer_size);
  fprintf(stdout,"table created: %d\n",1 << hash_key_bits);


  long count_files = 0;
  long long bad_reads =0;
  long long total_kmers =0;
  long long total_length = 0; //total sequence length


  //Go through all the files, loading data into the graph


  while (!feof(fp_fnames)){

    fscanf(fp_fnames, "%s\n", filename);
    
    count_files++;

    total_length += load_fasta_data_from_filename_into_graph(filename, &total_kmers, &bad_reads, MAX_READ_LENGTH, db_graph );

    fprintf(stdout,"\n\n************\n\n%ld file name:%s total_kmers: %lld, bad_reads: %lld, total seq:%lld\n\n",count_files,filename,total_kmers, bad_reads, total_length);

    //add data on current memory usage:
    FILE* fp_proc=fopen("/proc/self/status", "r");

    char line[500];
    while ( fgets(line,500,fp_proc) !=NULL)
      {
	if (line[0]=='V')
	  {
	  fprintf(stdout, "%s",line);
	  }
      }
    
  }

  printf("print supernodes\n");
  hash_table_traverse(&print_supernode,db_graph);

  return 1;
}
