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
  genome_complexity.c
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fnmatch.h>
#include <limits.h>

// cortex_var headers
#include "db_complex_genotyping.h"
#include "genome_complexity.h"
#include "db_variants.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "file_reader.h"
#include "maths.h"


boolean acgt(char c){
  boolean ret;
  if (c  != 'A' && c != 'a' && 
      c != 'C' && c != 'c' && 
      c != 'G' && c != 'g' && 
      c != 'T' && c != 't')
    {
      ret = false;
    }	
  else{
    ret =  true;
  }
  
  return ret;
}

boolean does_allele_lie_in_graph(dBNode** array_nodes, int len, int colour_cleaned_genome)
{
  int i;
  boolean ret = true;
  for (i=0; (ret==true) && (i<len); i++)
    {
      if (db_node_is_this_node_in_this_person_or_populations_graph(array_nodes[i], 
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
boolean allele_is_clean(dBNode** array_nodes,
			int len, int colour_cleaned_genome, int kmer)
{
  if ( (colour_cleaned_genome<0) || (colour_cleaned_genome>=NUMBER_OF_COLOURS) )
    {
      die("Cortex has been compiled for %d colours, but you are passing colour numnber %d into allele_is_clean\n", NUMBER_OF_COLOURS, colour_cleaned_genome);
    }
  int i;
  boolean all_nodes_look_ok=true;
  for (i=0; (i<len) && (i<=kmer) ; i++)//interior nodes
    {
      if (db_node_is_this_node_in_this_person_or_populations_graph(array_nodes[i],
								   colour_cleaned_genome)==true)
	{
	  Nucleotide nuc;
	  if (db_node_has_precisely_one_edge(array_nodes[i], forward, &nuc, colour_cleaned_genome)==false)
	    {
	      all_nodes_look_ok=false;
	    }
	}
    }

  return all_nodes_look_ok;
}

//given a fasta file, get one read at a time, and take a chunk 2*k+1 long from it,
//that becomes our notional ref-allele. If the ref-allele does not lie in the graph, ignore this read.
//then randomly mutate the central base, and that makes the alt allele
//If both ref and alt alleles are clean supernodes, increment total_errors_forming_clean_bubbles
void count_reads_where_snp_makes_clean_bubble(dBGraph* db_graph, char* fasta, boolean allow_reads_shorter_than_2k_plus_one, 
					      int colour_cleaned_genome,
					      int* total_errors_tested, int* total_errors_forming_clean_bubbles,
					      int (*file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
					      dBNode** array_nodes, Orientation* array_or, //assume these are length max_read_length+k+1 - plenty of space
					      Sequence* seq, KmerSlidingWindow* kmer_window, int max_read_length)				  
{
  FILE* fp = fopen(fasta, "r");
  if (fp == NULL){
    die("estimate_genome_complexity cannot open file:%s\n",fasta);
  }

  srand ( time(NULL) );

  int num_kmers_read=1;
  while ( (num_kmers_read>0) && (*total_errors_tested < MAX_NUM_READS_USED_FOR_ESTIMATE) )
    {
      boolean f_entry=true;
      num_kmers_read = align_next_read_to_graph_and_return_node_array(fp, max_read_length, array_nodes, array_or, 
								      false, &f_entry, file_reader, seq, kmer_window, db_graph, -1);

      if (!f_entry)
	{
	  die("One of these SNP reads is longer than the specified max read length\n");
	}
      if(num_kmers_read > 0 &&
         strlen(seq->seq) > (unsigned)db_graph->kmer_size/2 &&
         fnmatch("N", seq->seq, 0) != 0)
	{
    //not end of file & reasonable length read

	  boolean no_N_in_first_2k_plus_1=true;
	  int j;
	  for (j=0; (j<2*db_graph->kmer_size +1) && (j<(signed)strlen(seq->seq)); j++)
	    {
	      if (acgt((seq->seq)[j])==false)
		{
		  no_N_in_first_2k_plus_1=false; 
		}
	    }



	  if ((num_kmers_read<db_graph->kmer_size) && (allow_reads_shorter_than_2k_plus_one==false) )
	    {
	      //ignore this read - not long enough for a full SNP branch at this k
	    }
	  else if (no_N_in_first_2k_plus_1==false)
	    {
	      //ignore
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
		      
		      if (allele_is_clean(array_nodes, len, colour_cleaned_genome, db_graph->kmer_size)==false)
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
			  //int nkmers =
          get_single_kmer_sliding_window_from_sequence(alt,strlen(alt),
										    db_graph->kmer_size, 
										    kmer_window, db_graph);
			  //work through the sliding window and put nodes into the array you pass in. 
			  //Note this may find NULL nodes if the kmer is not in the graph
			  //also note seq is not used
			  load_kmers_from_sliding_window_into_array(kmer_window, 
								    db_graph, array_nodes, array_or, 
								    2*db_graph->kmer_size+1, false, -1);
			  if (allele_is_clean(array_nodes, len, colour_cleaned_genome, db_graph->kmer_size)==true)
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


/*
//assume graph contains cleaned genome, plus empty colour for working with
//will only use reads where entire read (ar at least the 2k bases around the error) lies in the pre-existing graph
double estimate_genome_complexity(dBGraph* db_graph, char* filename_fastaq,
				  boolean allow_reads_shorter_than_2k_plus_one, 
				  int colour_cleaned_genome, 
				  int max_read_length, FileFormat format,
				  int fastq_ascii_offset
				  )
{

  // ===========================================================================
  //   initialise stuff for reading of fasta;
  //     - since we support reading a LIST of fasta/q
  //       we alloc/init all this stuff then pass it
  //       down into the function that reads a single fasta
  // ===========================================================================
  

  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window in estimate_genome_complexity. Exit.\n");
    }
  

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer in estimate_genome_complexity. Exit.\n");
    }
  kmer_window->nkmers=0;
  
  
  //end of intialisation 
	  
	  
  //create file readers
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
      //die("new_entry must be true in hsi test function");
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }

  int file_reader_fastq(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      die("new_entry has to be true for fastq");
    }

    return read_sequence_from_fastq(fp,seq,max_read_length,fastq_ascii_offset);
  }

  

  if (max_read_length>MAX_READLEN_FOR_GEN_COMPLEXITY)
    {
      die("estime_genome_complexity is set up to handle short reads, and you have used max_read_lengh %d - not permitted\n", max_read_length);
    }
  dBNode* array_nodes[max_read_length+db_graph->kmer_size+1];
  Orientation array_or[max_read_length+db_graph->kmer_size+1];


  // ***********
  // end of init/malloc
  // ***********





  int total_errors_tested=0;//num of "alt alleles" tested 
  int total_errors_form_clean_bubbles=0; //number of these where both ref and alt allele form clean supernodes,
  if (format==FASTA)
    {
      count_reads_where_snp_makes_clean_bubble(db_graph, filename_fastaq, allow_reads_shorter_than_2k_plus_one, 
					       colour_cleaned_genome, 
					       &total_errors_tested, &total_errors_form_clean_bubbles,
					       file_reader_fasta, array_nodes, array_or, 
					       seq, kmer_window, max_read_length);
    }
  else if (format==FASTQ)
    {
      count_reads_where_snp_makes_clean_bubble(db_graph, filename_fastaq, allow_reads_shorter_than_2k_plus_one, 
					       colour_cleaned_genome,
					       &total_errors_tested, &total_errors_form_clean_bubbles,
					       file_reader_fastq, array_nodes, array_or, 
					       seq, kmer_window, max_read_length);
      
    }
 
  if (total_errors_tested==0)
    {
      printf("Unable to estimate genome complexity, returning zero. There are various possible reasons\n");
      printf("a) none of the reads in the fastaq you listed lay within the graph. That only happens if you specify fastaq/binary that do not correspond to each other, or possibly if you have done excessively brutal cleaning on the graph\n");
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

*/

// removed 'int colour_cleaned_genome' -- unused parameter
double estimate_genome_complexity(dBGraph* db_graph, char* fastaq,
                                  boolean allow_reads_shorter_than_2k_plus_one, 
                                  int working_colour, int max_read_length,
                                  FileFormat format, int fastq_ascii_offset,
                                  int* number_reads_used_in_estimate)
{



  //***********************************************
  //   initialise stuff for reading of fastaq;
  //***********************************************
  
  *number_reads_used_in_estimate=0;

  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window in estimate_genome_complexity. Exit.\n");
    }
  

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size+2));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer in estimate_genome_complexity. Exit.\n");
    }
  kmer_window->nkmers=0;

  //int max_windows = max_read_length/(db_graph->kmer_size+1);
  //number of possible kmers in a 'perfect' read
  //int max_kmers   = max_read_length-db_graph->kmer_size+1;

  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    die("Out of memory trying to allocate a KmerArraySet");
  }  
  windows->window=kmer_window;
  windows->nwindows=1;
  
  
  //end of intialisation 
	  
	  
  //create file readers
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
      //die("new_entry must be true in hsi test function");
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }

  int file_reader_fastq(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      die("new_entry has to be true for fastq");
    }

    return read_sequence_from_fastq(fp,seq,max_read_length,fastq_ascii_offset);
  }

  int (*file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry);
  if (format==FASTA)
    {
      file_reader=file_reader_fasta;
    }
  else if (format==FASTQ)
    {
      file_reader=file_reader_fastq;
    }
  else
    {
      die("Bad file format - in gen comp.\n");
    }
  

  if (max_read_length>MAX_READLEN_FOR_GEN_COMPLEXITY)
    {
      die(
"estime_genome_complexity is set up to handle relatively short reads, and you \n"
"have used max_read_lengh %d - not permitted.  Remove any reads longer than %d\n"
"and retry. Sorry - this is not very graceful\n",
          max_read_length, MAX_READLEN_FOR_GEN_COMPLEXITY);
    }

  


  dBNode** array_nodes = (dBNode**) malloc(sizeof(dBNode*)*(max_read_length+db_graph->kmer_size+1));
  Orientation* array_or= (Orientation*) malloc(sizeof(Orientation)*(max_read_length+db_graph->kmer_size+1));  
  dBNode** array_nodes_mut = (dBNode**) malloc(sizeof(dBNode*)*(max_read_length+db_graph->kmer_size+1));
  Orientation* array_or_mut= (Orientation*) malloc(sizeof(Orientation)*(max_read_length+db_graph->kmer_size+1));  

  dBNode** p1_nodes  = (dBNode**) malloc(sizeof(dBNode*)*(max_read_length+db_graph->kmer_size+1));
  Orientation* p1_or = (Orientation*) malloc(sizeof(Orientation)*(max_read_length+db_graph->kmer_size+1));  
  Nucleotide* p1_lab = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_read_length+db_graph->kmer_size+1));
  char* p1_str       = (char*) malloc(sizeof(char)*(max_read_length+db_graph->kmer_size+1) );
  dBNode** p2_nodes  = (dBNode**) malloc(sizeof(dBNode*)*(max_read_length+db_graph->kmer_size+1));
  Orientation* p2_or = (Orientation*) malloc(sizeof(Orientation)*(max_read_length+db_graph->kmer_size+1));  
  Nucleotide* p2_lab = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_read_length+db_graph->kmer_size+1));
  char* p2_str       = (char*) malloc(sizeof(char)*(max_read_length+db_graph->kmer_size+1) );

  if ( (array_nodes==NULL) || (array_or==NULL) || (array_nodes_mut==NULL) || (array_or_mut==NULL) 
       || (p1_nodes==NULL) || (p1_or==NULL) || (p1_lab==NULL) || (p1_str==NULL)
       || (p2_nodes==NULL) || (p2_or==NULL) || (p2_lab==NULL) || (p2_str==NULL) )
    {
      die("Cannot malloc node arrays for estimating genome complexity - you must be very short of memor\n");
    }
       
  /*
  dBNode* array_nodes[max_read_length+db_graph->kmer_size+1];
  Orientation array_or[max_read_length+db_graph->kmer_size+1];
  dBNode* array_nodes_mut[max_read_length+db_graph->kmer_size+1];
  Orientation array_or_mut[max_read_length+db_graph->kmer_size+1];
  dBNode* p1_nodes[max_read_length+db_graph->kmer_size+1];
  Orientation p1_or[max_read_length+db_graph->kmer_size+1];
  Nucleotide p1_lab[max_read_length+db_graph->kmer_size+1];
  char p1_str[max_read_length+db_graph->kmer_size+1];
  dBNode* p2_nodes[max_read_length+db_graph->kmer_size+1];
  Orientation p2_or[max_read_length+db_graph->kmer_size+1];
  Nucleotide p2_lab[max_read_length+db_graph->kmer_size+1];
  char p2_str[max_read_length+db_graph->kmer_size+1];
  */
  
  //will ignore contents of this - needed for API - get read-len distribution (after filters):
  long long* readlen_distrib=(long long*) malloc(sizeof(long long) * (max_read_length+1));
  long long** readlen_distrib_ptrs = (long long**) malloc(sizeof(long long*) * (max_read_length+1));
  
  if ( (readlen_distrib==NULL) || (readlen_distrib_ptrs==NULL) )
    {
      die("Unable to malloc array to hold readlen distirbution!Exit.\n");
    }
  int i;
  for (i=0; i<=max_read_length; i++)
    {
      readlen_distrib[i]=0;
      readlen_distrib_ptrs[i]=&readlen_distrib[i];
    }
  
  // ***********
  // end of init/malloc
  // ***********

  //  if (genome_graph_built_from_seq_data==true)
  //  {
  //       db_graph_smooth_bubbles_for_specific_person_or_pop(6,db_graph->kmer_size+1, db_graph, colour_cleaned_genome);
      //  }
  

  //load the fasta into the working colour
  /*  
  if (format==FASTA)
    {
      load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(fastaq, &bases_read, &bases_loaded, &bad_reads, &dup_reads, max_read_length,
									 remove_duplicates_single_endedly, break_homopolymers, homopol, 
									 db_graph, working_colour);
    }
  else
    {
      //get read-len distribution (after filters):
      long long* readlen_distrib=(long long*) malloc(sizeof(long long) * (max_read_length+1));
      long long** readlen_distrib_ptrs = (long long**) malloc(sizeof(long long*) * (max_read_length+1));

      if ( (readlen_distrib==NULL) || (readlen_distrib_ptrs==NULL) )
	{
	  die("Unable to malloc array to hold readlen distirbution!Exit.\n");
	}
      int i;
      for (i=0; i<=max_read_length; i++)
	{
	  readlen_distrib[i]=0;
	  readlen_distrib_ptrs[i]=&readlen_distrib[i];
	}


      
      load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(fastaq, &bases_read, &bases_loaded,
									 readlen_distrib_ptrs,  &bad_reads, 0, &dup_reads, 
									 max_read_length,
									 remove_duplicates_single_endedly, break_homopolymers,
									 homopol,fastq_ascii_offset, 
									 db_graph, working_colour);
      free(readlen_distrib);
      free(readlen_distrib_ptrs);
    }
  */

  int total_errors_tested=0;//num of "alt alleles" tested 
  int total_errors_form_clean_bubbles=0; //number of these where both ref and alt allele form clean supernodes,


  FILE* fp = fopen(fastaq, "r");
  if (fp == NULL){
    die("estimate_genome_complexity  cannot open file:%s\n",fastaq);
  }

  srand ( time(NULL) );

  int num_kmers_read=1;
  while ( (num_kmers_read>0) && (total_errors_tested < MAX_NUM_READS_USED_FOR_ESTIMATE) )
    {
      
      boolean f_entry=true;
      num_kmers_read = align_next_read_to_graph_and_return_node_array(fp, max_read_length, array_nodes, array_or, 
								      false, &f_entry, file_reader, seq, kmer_window, db_graph, -1);

      if (!f_entry)
	{
	  die("One of these reads is longer than the specified max read length\n");
	}
      if ((num_kmers_read>0)&&(strlen(seq->seq)>1+(unsigned)db_graph->kmer_size /2) ) 
	//not end of file & reasonable length read 
	{
	  boolean no_N_in_first_2k_plus_1=true;
	  boolean quality_above_10=true;
	  int j;
	  for (j=0; (j<2*db_graph->kmer_size +1) && (j<(signed)strlen(seq->seq)); j++)
	    {
	      if (acgt((seq->seq)[j])==false)
		{
		  no_N_in_first_2k_plus_1=false; 
		}
	      //if (seq->qual[j]<10)
	      //	{
	      //	  quality_above_10=false;
	      //	}
	    }
	  


	  if ((num_kmers_read<db_graph->kmer_size+1) && (allow_reads_shorter_than_2k_plus_one==false) )
	    {
	      //ignore this read - not long enough for a full SNP branch at this k plus a node before that to be 
	      //the junction at the start of both branches
	      //printf("Ignore this read %s\n", seq->seq);
	    }
	  else if (no_N_in_first_2k_plus_1==false)
	    {
	      //ignore
	      //printf("Ignore this read %s\n", seq->seq);

	    }
	  else if (quality_above_10==false)
	    {
	    }
	  else
	    {
	      //printf("this is good: %s\n", seq->seq);
	      //make a mutated allele
	      char alt[2*db_graph->kmer_size+2];
	      int j;
	      for (j=0; (j<2*db_graph->kmer_size+1) && (j<(signed)strlen(seq->seq)); j++)
		{
		  alt[j]= (seq->seq)[j];
		}
	      alt[j]='\0';
	      //strcpy(alt, seq->seq);
	      int char_to_modify=db_graph->kmer_size; //ie leave one base to be the same for both branches, then k-th base after that
	               
	      modify_character(alt, char_to_modify, rand() % 3 );
	      
	      //turn it into a sliding window 
	      //int nkmers =
          get_single_kmer_sliding_window_from_sequence(alt,strlen(alt),
									db_graph->kmer_size, 
									kmer_window, db_graph);

	      boolean prev_full=true;
	      long long b_loaded=0;

	      load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(
          windows, &prev_full, &b_loaded, false, db_graph,
          working_colour, readlen_distrib_ptrs);

	      //work through the sliding window and put nodes into the array you pass in. 
	      //Note this may find NULL nodes if the kmer is not in the graph
	      //also note seq is not used
	      load_kmers_from_sliding_window_into_array(kmer_window, 
							db_graph, array_nodes_mut, array_or_mut, 
							2*db_graph->kmer_size+1, false, -1);
	      double av1, av2;
	      Covg min1,min2,max1,max2;
	      int len1=0;
	      int len2=0;
	      boolean forms_clean_bubble=db_graph_detect_bubble_in_subgraph_defined_by_func_of_colours(&db_node_condition_always_true,
												       array_nodes_mut[0],
												       array_or_mut[0], db_graph->kmer_size+1,
												       &db_node_action_do_nothing,
												       &len1, p1_nodes, p1_or, p1_lab, p1_str,
												       &av1, &min1, &max1,
												       &len2,p2_nodes, p2_or, p2_lab, p2_str,
												       &av2, &min2, &max2,
												       db_graph, &element_get_colour_union_of_all_colours,
												       &element_get_covg_union_of_all_covgs, false, NULL);

	      total_errors_tested++;
	      *number_reads_used_in_estimate = *number_reads_used_in_estimate +1;
	      if (forms_clean_bubble==true) 
		{
		  total_errors_form_clean_bubbles++;
		  //printf("Good read %s, this error %s is clean\n", seq->seq, alt);
		}
	      else
		{
		  //printf("Bad read %s, this error %s is not clean\n", seq->seq, alt);
		}
	    } 
	}
    }

  if (total_errors_tested==0)
    {
      printf("Unable to estimate genome complexity, returning zero. There are various possible reasons\n");
      printf("a) none of the reads in the fastaq you listed lay within the graph. That only happens if you specify fastaq/binary that do not correspond to each other, or possibly if you have done excessively brutal cleaning on the graph\n");
      if (allow_reads_shorter_than_2k_plus_one==false)
	{
	  printf("b) Your fastaq has reads <2*kmer length, so you can't put an error in the middle of a read and work out a full bubble branch. See Manual and our paper\n");
	}
      return (double) 0;
    }

  //cleanup
  free_sequence(&seq);
  free(windows->window->kmer);
  free(windows->window);
  free(windows);
  free(readlen_distrib);
  free(readlen_distrib_ptrs);
  free(array_nodes);
  free(array_or);
  free(array_nodes_mut);
  free(array_or_mut);
  free(p1_nodes);
  free(p1_or);
  free(p1_lab);
  free(p1_str);
  free(p2_nodes);
  free(p2_or);
  free(p2_lab);
  free(p2_str);
  
  return (double)total_errors_form_clean_bubbles / (double) total_errors_tested;
}
