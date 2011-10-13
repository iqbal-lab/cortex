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

#include <db_complex_genotyping.h>
#include <fnmatch.h>
#include <genome_complexity.h>
#include <db_variants.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <global.h>
#include <file_reader.h>
#include <maths.h>
#include <limits.h>


boolean does_allele_lie_in_graph(dBNode** array_nodes, int len, int colour_cleaned_genome)
{
  int i;
  boolean ret = true;
  for (i=0; (ret==true) && (i<len); i++)
    {
      if (db_node_is_this_node_in_this_person_or_populations_graph(array_nodes[i], 
								   individual_edge_array, 
								   colour_cleaned_genome)==false)
	{
	  ret=false;
	}
    }
  return ret;
}

//I don't care if the nodes exist in the hash, nor if it is in this colour. 
//but if it is in the graph and in the colour, none of the interior nodes are 
//allowed to have in/out degree>1
boolean allele_is_clean(dBNode** array_nodes,Orientation* array_or, 
			int len, int colour_cleaned_genome, int kmer)
{
  if ( (colour_cleaned_genome<0) || (colour_cleaned_genome>=NUMBER_OF_COLOURS) )
    {
      printf("Cortex has been compiled for %d colours, but you are passing colour numnber %d into allele_is_clean\n", NUMBER_OF_COLOURS, colour_cleaned_genome);
      exit(1);
    }
  int i;
  boolean all_interior_nodes_look_ok=true;
  for (i=1; (i<len) && (i<kmer) ; i++)//interior nodes
    {
      if (db_node_is_this_node_in_this_person_or_populations_graph(array_nodes[i],
								   individual_edge_array,
								   colour_cleaned_genome)==true)
	{
	  Nucleotide nuc;
	  if (db_node_has_precisely_one_edge(array_nodes[i], forward, &nuc, individual_edge_array, colour_cleaned_genome)==false)
	    {
	      all_interior_nodes_look_ok=false;
	    }
	}
    }

  return all_interior_nodes_look_ok;
}

//given a fasta file, get one read at a time, and take a chunk 2*k+1 long from it,
//that becomes our notional ref-allele. If the ref-allele does not lie in the graph, ignore this read.
//then randomly mutate the central base, and that makes the alt allele
//If both ref and alt alleles are clean supernodes, increment total_errors_forming_clean_bubbles
void get_clean_and_unclean_counts(dBGraph* db_graph, char* fasta, boolean allow_reads_shorter_than_2k_plus_one, 
				  int colour_cleaned_genome,
				  int* total_errors_tested, int* total_errors_forming_clean_bubbles,
				  int (*file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
				  dBNode** array_nodes, Orientation* array_or, //assume these are length max_read_length+k+1 - plenty of space
				  Sequence* seq, KmerSlidingWindow* kmer_window, int max_read_length)				  
{
  FILE* fp = fopen(fasta, "r");
  if (fp == NULL){
    printf("estimate_genome_complexity  cannot open file:%s\n",fasta);
    exit(1);
  }

  srand ( time(NULL) );

  int num_kmers_read=1;
  while ( (num_kmers_read>0) && (*total_errors_tested < MAX_NUM_READS_USED_FOR_ESTIMATE) )
    {
      num_kmers_read = align_next_read_to_graph_and_return_node_array(fp, max_read_length, array_nodes, array_or, 
								      false, file_reader, seq, kmer_window, db_graph, -1);

      
      if ((num_kmers_read>0)&&(strlen(seq->seq)>db_graph->kmer_size /2)&& (fnmatch("N", seq->seq, 0)!=0) ) 
	//not end of file & reasonable length read & no N's in read
	{
	  if ((num_kmers_read<db_graph->kmer_size) && (allow_reads_shorter_than_2k_plus_one==false) )
	    {
	      //ignore this read - not long enough for a full SNP branch at this k
	    }
	  else
	    {
	      if (does_allele_lie_in_graph(array_nodes, num_kmers_read, colour_cleaned_genome)==true)
		{
		  //we want to look at the first k kmers ideally, but we may have less than that number of kmers (short read length compared with k, or Ns for example)
		  int len = min_of_ints(num_kmers_read, db_graph->kmer_size);
		  
		  if (len>db_graph->kmer_size/2)//a bit arbitrary
		    {
		      //this is a testable site/use-able read.
		      *total_errors_tested = *total_errors_tested +1;
		      
		      if (allele_is_clean(array_nodes, array_or, len, colour_cleaned_genome, db_graph->kmer_size)==false)
			{
			  // error does not make a clean bubble
			}
		      else
			{
			  //make a mutated allele
			  char alt[2*db_graph->kmer_size+1];
			  int j;
			  for (j=0; j<2*db_graph->kmer_size+1; j++)
			    {
			      alt[j]='\0';
			    }
			  strcpy(alt, seq->seq);
			  int char_to_modify=db_graph->kmer_size -1; //ie k-th base
			  modify_character(alt, char_to_modify, rand() % 3 );
			  
			  //turn it into a sliding window 
			  int nkmers = get_single_kmer_sliding_window_from_sequence(alt,strlen(alt),
										    db_graph->kmer_size, 
										    kmer_window, db_graph);
			  //work through the sliding window and put nodes into the array you pass in. 
			  //Note this may find NULL nodes if the kmer is not in the graph
			  //also note seq is not used
			  load_kmers_from_sliding_window_into_array(kmer_window, seq, 
								    db_graph, array_nodes, array_or, 
								    len, false, -1);
			  if (allele_is_clean(array_nodes, array_or, len, colour_cleaned_genome, db_graph->kmer_size)==true)
			    {
			      *total_errors_forming_clean_bubbles = *total_errors_forming_clean_bubbles+1; 
			    }

			}
		    }
		}
	    }
	  
	}

    }
  

  
}



//assume graph contains cleaned genome, plus empty colour for working with
//will only use reads where entire read (ar at least the 2k bases around the error) lies in the pre-existing graph
double estimate_genome_complexity(dBGraph* db_graph, char* filename_list_fastaq,
				  boolean allow_reads_shorter_than_2k_plus_one, 
				  int colour_cleaned_genome, 
				  int max_read_length, FileFormat format,
				  int fastq_ascii_offset
				  )
{

  //***********************************************
  //   initialise stuff for reading of fasta;
  //     - since we support reading a LIST of fasta/q
  //       we alloc/init all this stuff then pass it
  //       down into the function that reads a single fasta
  //***********************************************
  

  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in estimate_genome_complexity. Exit.\n");
      exit(1);
    }
  

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in estimate_genome_complexity. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  
  
  //end of intialisation 
	  
	  
  //create file readers
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
      //printf("new_entry must be true in hsi test function");
      //exit(1);
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }

  int file_reader_fastq(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      puts("new_entry has to be true for fastq\n");
      exit(1);
    }

    return read_sequence_from_fastq(fp,seq,max_read_length,fastq_ascii_offset);
  }

  

  if (max_read_length>MAX_READLEN_FOR_GEN_COMPLEXITY)
    {
      printf("estime_genome_complexity is set up to handle short reads, and you have used max_read_lengh %d - not permitted\n", max_read_length);
      exit(1);
    }
  dBNode* array_nodes[max_read_length+db_graph->kmer_size+1];
  Orientation array_or[max_read_length+db_graph->kmer_size+1];


  // ***********
  // end of init/malloc
  // ***********





  FILE* fp = fopen(filename_list_fastaq, "r");
  if (fp == NULL){
    printf("estimate_genome_complexity  cannot open file:%s\n",filename_list_fastaq);
    exit(1);
  }

  int total_errors_tested=0;//num of "alt alleles" tested 
  int total_errors_form_clean_bubbles=0; //number of these where both ref and alt allele form clean supernodes,
  char fastaq[MAX_FILENAME_LENGTH+1];


  while(fgets(fastaq,MAX_FILENAME_LENGTH, fp) !=NULL)
    {

      //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(fastaq, '\n')) != NULL)
	{
	  *p = '\0';
	}

      if (format==FASTA)
	{
	  get_clean_and_unclean_counts(db_graph, fastaq, allow_reads_shorter_than_2k_plus_one, 
				       colour_cleaned_genome, 
				       &total_errors_tested, &total_errors_form_clean_bubbles,
				       file_reader_fasta, array_nodes, array_or, 
				       seq, kmer_window, max_read_length);
	}
      else if (format==FASTQ)
	{
	  get_clean_and_unclean_counts(db_graph, fastaq, allow_reads_shorter_than_2k_plus_one, 
				       colour_cleaned_genome,
				       &total_errors_tested, &total_errors_form_clean_bubbles,
				       file_reader_fastq, array_nodes, array_or, 
				       seq, kmer_window, max_read_length);
	  
	}
 
    }
  close(fp);

  if (total_errors_tested==0)
    {
      printf("Unable to estimate genome complexity, returning zero. There are various possible reasons\n");
      printf("a) none of the reads in the fastaq you listed lay within the graph. That only happens if you specify fastaq/binary that do not sorrespond to each other, or possibly if you have done excessively brutal cleaning on the graph\n");
      if (allow_reads_shorter_than_2k_plus_one==false)
	{
	  printf("b) Your fastaq has reads <2*kmer length, so you can't put an error in the middle of a read and work out a full bubble branch. See Manual and our paper\n");
	}
      return (double) 0;
    }

  //cleanup
  free_sequence(&seq);
  free(kmer_window->kmer);
  free(kmer_window);

  return (double)total_errors_form_clean_bubbles / (double) total_errors_tested;
}
