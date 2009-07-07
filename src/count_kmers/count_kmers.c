#include <seq.h>
#include <binary_kmer.h>
#include <open_hash/hash_table.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){
  Sequence * seq;
  FILE *fp_fnames,*fp_file;
 
  char filename[100];
  int kmer_size, hash_key_bits;

  HashTable * hash_table = NULL;
   
  long count=0;
  FILE *fout;
  
  void print_element_count(Element * e){
  
    char filename [50];
    if (count % 100000000 == 0){
      int index = count / 100000000;
      
      if (count !=0){
        fclose(fout);
      }
      
      sprintf(filename,"out_kmers_%i_%i",kmer_size,index);
      
      fprintf(stderr,"opening file %s\n",filename);
      fout = fopen(filename,"w");
    }
    
    count++;
    
    element_print(fout,e,hash_table->kmer_size,NULL);
    
  }




  kmer_size        = atoi(argv[2]);
  hash_key_bits    = atoi(argv[3]); //2 ** hash_key_bits buckets

  fprintf(stderr,"Kmer size:%d hash_table_size (bits):%d\n",kmer_size,hash_key_bits);

  hash_table = hash_table_new(hash_key_bits,kmer_size);

  fprintf(stderr,"hash table created with this many buckets: %i\n",1<<hash_key_bits);

  //open file of file names
  fp_fnames= fopen(argv[1], "r");

  int count_kmer   = 0;
  int count_file   = 0;
  long long total_length = 0; //total sequence length

  while (!feof(fp_fnames)){

    int count_bad_reads = 0;
    fscanf(fp_fnames, "%s\n", filename);
    
    fp_file = fopen(filename, "r");
    if (fp_file == NULL){
      printf("cannot open file:%s\n",filename);
      exit(1);
    }
    
    int seq_length = 0;
    count_file++;

    while ((seq = read_sequence_from_fasta(fp_file))){
      KmerArray *kmers;
      int i;

      seq_length += seq->length;

      kmers = get_binary_kmers_from_sequence(seq->seq,seq->length,kmer_size);
      free_sequence(&seq); 

      if (kmers == NULL){
	count_bad_reads++;
      }
      else{
	
	for(i=0;i<kmers->nkmers;i++){
	  hash_table_apply_or_insert(element_get_key(kmers->bin_kmers[i],hash_table->kmer_size),&element_increment_count,hash_table);
	  
      }
      
	binary_kmer_free_kmers(&kmers);
      }
    }  
     
    total_length+=seq_length;
    fprintf(stderr,"%i file name:%s kmers:%i bad_reads:%i seq:%i total seq:%qd\n",count_file,filename,count_kmer,count_bad_reads,seq_length, total_length);
    
  }


  fprintf(stderr,"printing table\n");
  hash_table_traverse(&print_element_count,hash_table);

  return 0;
}
  
