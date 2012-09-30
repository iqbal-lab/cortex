/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo    
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  count_kmers.c - not sure if this file is still used -- doesn't appear to compile
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

// cortex_var headers
#include "seq.h"
#include "binary_kmer.h"
#include "open_hash/hash_table.h"

int main(int argc, char **argv){
  Sequence * seq;
  FILE *fp_fnames,*fp_file;
 
  char filename[100];
  int max_read_length = 2000;
  int fastq; //if 0 entry is fasta otherwise is quality cut-off
  int kmer_size, hash_key_bits;
  int bucket_size;
  HashTable * hash_table = NULL;
  long long entry_length;

  boolean full_entry = true;
  boolean prev_full_entry = true;
  Element * prev_entry =  NULL;

  long count=0;
  FILE *fout;
  

  //useful functions
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


  int fastq_file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      die("new_entry has to be true for fastq");
    }

    return read_sequence_from_fastq(fp,seq,max_read_length);
  }
  

  int fasta_file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = kmer_size;
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }


  
  //action starts here
  kmer_size        = atoi(argv[2]);
  hash_key_bits    = atoi(argv[3]); //2 ** hash_key_bits buckets
  bucket_size      = atoi(argv[4]);
  fastq            = atoi(argv[5]);

  if (kmer_size>31){
    die("kmer size [%i] exceeds 31");
  }

  fprintf(stderr,"Input file of filenames: %s\n",argv[1]);
  fprintf(stderr,"Kmer size: %d hash_table_size (%d bits): %d - bucket size: %d - total size: %qd\n",kmer_size,hash_key_bits,1 << hash_key_bits, bucket_size, ((long long) 1<<hash_key_bits)*bucket_size);


  hash_table = hash_table_new(hash_key_bits,bucket_size, 10,kmer_size); //10 is the number of rehashing attemps 
  fprintf(stderr,"table created: %d\n",1 << hash_key_bits); 

  if (fastq>0){
    fprintf(stderr,"fastq files: quality cut-off: %i\n",fastq);
  }
  else{
    fprintf(stderr,"fasta files\n");
  }

  //open file of file names
  fp_fnames= fopen(argv[1], "r");


  int count_file   = 0;
  long long total_length = 0; //total sequence length


  //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  
 


  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call 
  //-- with the view to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    die("Out of memory trying to allocate a KmerArraySet");
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);
  

  //for each file 
  while (!feof(fp_fnames)){

    long long count_bad_reads = 0;
    fscanf(fp_fnames, "%s\n", filename);
    
    fp_file = fopen(filename, "r");
    if (fp_file == NULL){
      die("cannot open file:%s",filename);
    }
    
    long long seq_length = 0;
    count_file++;
    
    while ((entry_length = fastq == 0 ? 
			   fasta_file_reader(fp_file,seq,max_read_length,full_entry,&full_entry) :
			   fastq_file_reader(fp_file,seq,max_read_length,full_entry,&full_entry))){

      int i,j;

      seq_length += (long long) (entry_length - (prev_full_entry==false ? kmer_size : 0));
   
      int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,
						     entry_length,fastq,kmer_size,windows,max_windows, max_kmers);

      
      if (nkmers == 0) {
	count_bad_reads++;
      }
      else {

	 for(i=0;i<windows->nwindows;i++){ //for each window
	   KmerSlidingWindow * current_window = &(windows->window[i]);
	   
	   for(j=0;j<current_window->nkmers;j++){ //for each kmer in window
	     boolean found = false;

	     Element * current_entry = hash_table_find_or_insert(element_get_key(current_window->kmer[j],kmer_size),&found,hash_table);	
 
	     if (current_entry == NULL){
	       die("file_reader: problem - current kmer not found");
	     }
	     
	     if (! (i==0 && j==0 && prev_full_entry == false && current_entry == prev_entry)){ //otherwise is the same old last entry 
	       element_increment_count(current_entry);
	     }
	  
	   }
	 }
      }
      
      //for long fasta entries this pushes the last kmer back to the beginning of the sequence
      if (full_entry == false){
	shift_last_kmer_to_start_of_sequence(seq,entry_length,kmer_size);
      }
      
      prev_full_entry = full_entry;
      
    }
     
    total_length+=seq_length;
    
    fprintf(stderr,"\n%i kmers: %qd file name:%s bad reads: %qd seq:%qd total seq:%qd\n\n",count_file,hash_table_get_unique_kmers(hash_table),filename,count_bad_reads,seq_length, total_length);

    
  }


  fprintf(stderr,"printing table\n");
  hash_table_traverse(&print_element_count,hash_table);

  return 0;
}
  
