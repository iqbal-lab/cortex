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

#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <cmd_line.h>
#include <time.h>
#include <graph_info.h>
#include <db_differentiation.h>
#include <math.h>
#include <maths.h>
#include <seq_error_rate_estimation.h>
#include <string.h>

void estimate_seq_error_rate_from_snps_for_each_colour(char* colourlist_snp_alleles, GraphInfo* db_graph_info, dBGraph* db_graph, int ref_colour, long long genome_size,
						       long double default_seq_err_rate)
{

  int max_read_length = 2*(db_graph->kmer_size);

  // ****** malloc and initialisation //
  
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
      printf("Failed to malloc kmer sliding window in align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. Exit.\n");
      exit(1);
    }
  

  //  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size+1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  
  
	  
	  
  //create file reader
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

  dBNode** array_nodes = (dBNode**) malloc(sizeof(dBNode*) * (max_read_length+db_graph->kmer_size+1) );
  Orientation*  array_or = (Orientation*) malloc(sizeof(Orientation)*(max_read_length+db_graph->kmer_size+1) );
  if ( (array_nodes==NULL) || (array_or==NULL) )
    {
      printf("Unable to malloc arrays for sequencing eror rate estimation - is your server low on memory?\n");
      exit(1);
    }

  //end of initialisation 




  int i;

  printf("Open colourlist %s\n", colourlist_snp_alleles);
  FILE* fp = fopen(colourlist_snp_alleles, "r");
  if (fp==NULL)
    {
      printf("Cannot open %s\n", colourlist_snp_alleles);
      exit(1);
    }

  char filename[MAX_FILENAME_LENGTH];
  int colour=-1;
  while (feof(fp) ==0)
    {
      if (fgets(filename,MAX_FILENAME_LENGTH, fp) !=NULL)
        {
	  char* p;
	  if ((p = strchr(filename, '\n')) != NULL)
	    {
	      *p = '\0';
	    }
	  
	  colour++;
	  if (colour != ref_colour)
	    {
	      int num_snps_tested=0;
	      db_graph_info->seq_err[colour]=estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta(filename, db_graph, colour, 
												     seq, kmer_window, &file_reader_fasta,
												     array_nodes, array_or, &num_snps_tested,
												     max_read_length, default_seq_err_rate);

	    }

	}
    }
  fclose(fp);

  free_sequence(&seq);
  free(kmer_window->kmer);
  free(kmer_window);
  free(array_nodes);
  free(array_or);
}



//pass in a fasta of SNP alleles that we know from SNP-chip genotyping, this sample (colour) does not contain.
//format 
//   >allele with expected zero covg
//   GAAGAGAGAG
//   >allele which we think is homozygous in this sample
//   GAAGATAGAG
//  should be k-1 bases before and after the SNP base itself.
long double estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta(char* fasta, dBGraph* db_graph, int colour, 
									 Sequence* seq, KmerSlidingWindow* kmer_window,
									 int (*file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
									 dBNode** array_nodes, Orientation* array_or, int* num_snps_tested, int max_read_length,
									 long double default_seq_err)
{
  int num_snps=0;
  long long total_covg_on_true_alleles=0; //2nd read in each pair
  long long total_covg_on_error_alleles=0; //first read in each pair

  FILE* fptr = fopen(fasta, "r");
  if (fptr==NULL)
    {
      printf("Unable to open file Z%sZ. Abort.\n", fasta);
      exit(1);
    }
  
  int num_kmers=0;
  boolean full_entry=true;
  do
    {
      num_kmers = align_next_read_to_graph_and_return_node_array(fptr, max_read_length, array_nodes, array_or, false, &full_entry, file_reader,
								 seq, kmer_window, db_graph, -1);

      printf("Next SNP\n");
      boolean too_short = false;
      total_covg_on_error_alleles += count_reads_on_allele_in_specific_colour(array_nodes, num_kmers, colour, &too_short);
      printf("New total covg on bad alleles is %qd\n", total_covg_on_error_alleles);
      num_kmers = align_next_read_to_graph_and_return_node_array(fptr, max_read_length, array_nodes, array_or, false, &full_entry, file_reader,
								 seq, kmer_window, db_graph, -1);

      total_covg_on_true_alleles += count_reads_on_allele_in_specific_colour(array_nodes, num_kmers, colour, &too_short);
      printf("New total covg on true alleles is %qd\n", total_covg_on_true_alleles);

      num_snps++; 
    }
  while((num_kmers>0)||(!feof(fptr)) );

    
    if (num_snps==0)
      {
	printf("This file %s seems to contain no reads. Will forcibly set estimated seq error rate to %Lf and continue. But you should look into this\n",  fasta, default_seq_err);
	return (long double) default_seq_err;
      }
    else
      {
	*num_snps_tested = num_snps;
      }

    if (total_covg_on_true_alleles+total_covg_on_error_alleles==0)
      {
	printf("RET DEF\n");
	return default_seq_err;
      }
    else
      {
	return total_covg_on_error_alleles/(total_covg_on_true_alleles+total_covg_on_error_alleles);
      }

  fclose(fptr);

}
