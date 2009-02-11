#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>

int main(int argc, char **argv){

<<<<<<< local
  FILE *fp_fnames;
  char filename[100];
=======
  FILE *fp_fnames,*fp_file;
  char filename[300];
>>>>>>> other
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


  int count_file   = 0;
  long long total_length = 0; //total sequence length

  //Go through all the files, loading data into the graph

  int total_reallocs=0;

  while (!feof(fp_fnames)){

    //int count_bad_reads = 0;
    fscanf(fp_fnames, "%s\n", filename);
    
    
    int seq_length = 0;
    count_file++;

<<<<<<< local
    total_length += load_fasta_data_from_filename_into_graph(filename, db_graph);
=======
    int reallocs=0;
    total_length += load_fasta_data_into_graph_efficient(fp_file, db_graph, &reallocs);
    total_reallocs +=reallocs;
>>>>>>> other

    fprintf(stdout,"\n\n************\n\n%i file name:%s seq:%i total seq:%qd\n\n",count_file,filename,seq_length, total_length);

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

  printf("Total number of reallocs is %d", total_reallocs);
  printf("print supernodes\n");
  hash_table_traverse(&print_supernode,db_graph);

  return 1;
}
