/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo 
 * 
 * CORTEX project contacts:  
 *         M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 *         Z. Iqbal (zam@well.ox.ac.uk)
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
  seq_error_rate_estimation.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// cortex_var headers
#include "element.h"
#include "file_reader.h"
#include "cmd_line.h"
#include "graph_info.h"
#include "db_differentiation.h"
#include "maths.h"
#include "seq_error_rate_estimation.h"

void estimate_seq_error_rate_from_snps_for_each_colour(char* colourlist_snp_alleles, 
                                                       GraphInfo* db_graph_info, 
                                                       dBGraph* db_graph, 
                                                       int ref_colour, 
                                                       //long long genome_size,
                                                       long double default_seq_err_rate, 
                                                       char* output_file)
{

  int max_read_length = 2*(db_graph->kmer_size);

  // ****** malloc and initialisation //
  
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
    die("Failed to malloc kmer sliding window in align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. Exit.\n");
  }
  

  //  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size+1));
  if (kmer_window->kmer==NULL)
  {
    die("Failed to malloc kmer_window->kmer in align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. Exit.\n");
  }
  kmer_window->nkmers=0;
  
  


  //create file reader
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

  dBNode** array_nodes = (dBNode**) malloc(sizeof(dBNode*) * (max_read_length+db_graph->kmer_size+1) );
  Orientation*  array_or = (Orientation*) malloc(sizeof(Orientation)*(max_read_length+db_graph->kmer_size+1) );
  if ( (array_nodes==NULL) || (array_or==NULL) )
  {
    die("Unable to malloc arrays for sequencing eror rate estimation - is your server low on memory?\n");
  }

  FILE* fout=NULL;
  if (output_file != NULL)
  {
    fout = fopen(output_file, "w");
    if (fout==NULL)
    {
      die("Unable to open output file %s\n", output_file);
    }
  }
  //end of initialisation 



  // Get absolute path
  char absolute_path[PATH_MAX+1];
  char* filename_abs_path = realpath(colourlist_snp_alleles, absolute_path);

  if(filename_abs_path == NULL)
  {
    die("Cannot get absolute path to seq_error colours: %s\n", filename_abs_path);
  }

  // Get directory path
  StrBuf *dir = file_reader_get_strbuf_of_dir_path(filename_abs_path);

  printf("Open colourlist %s\n", colourlist_snp_alleles);

  FILE* fp = fopen(colourlist_snp_alleles, "r");
  if(fp == NULL)
  {
    die("Cannot open %s\n", colourlist_snp_alleles);
  }

  StrBuf *line = strbuf_new();

  int colour = -1;

  while(strbuf_reset_readline(line, fp))
  {
    strbuf_chomp(line);

    if(strbuf_len(line) > 0)
    {
      colour++;

      if(colour != ref_colour)
      {
        // Get paths relative to filelist dir
        if(strbuf_get_char(line, 0) != '/')
          strbuf_insert(line, 0, dir, 0, strbuf_len(dir));

        // Get absolute paths
        char* path_ptr = realpath(line->buff, absolute_path);

        if(path_ptr == NULL)
        {
          die("Cannot find sequence file for seq_error_estimation: %s\n",
            line->buff);
        }

        // unused
        //int num_snps_tested = 0;


        long double seq_err = estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta(
          path_ptr, db_graph, colour, seq, kmer_window, &file_reader_fasta,
          array_nodes, array_or, //&num_snps_tested,
          max_read_length, default_seq_err_rate, fout);

        db_graph_info->seq_err[colour] = seq_err;
      }
    }
  }

  // Cleanup
  strbuf_free(line);
  strbuf_free(dir);

  fclose(fp);

  if(output_file != NULL)
  {
    fclose(fout);
  }

  free_sequence(&seq);
  free(kmer_window->kmer);
  free(kmer_window);
  free(array_nodes);
  free(array_or);
}



// Pass in a fasta of SNP alleles that we know from SNP-chip genotyping, this
// sample (colour) does not contain. format 
//   >allele with expected zero covg
//   GAAGAGAGAG
//   >allele which we think is homozygous in this sample
//   GAAGATAGAG
//  should be k-1 bases before and after the SNP base itself.
// if you enter NULL for output filename, will not dump distribution of numbers
// of SNPs with 0,..100 covg on the allele we expect to have zero covg on
long double estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta(
  char* fasta, dBGraph* db_graph, int colour, 
  Sequence* seq, KmerSlidingWindow* kmer_window,
  int (*file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
  dBNode** array_nodes, Orientation* array_or, //int* num_snps_tested, //unused
  int max_read_length, long double default_seq_err, FILE* fout)
{

  int num_snps_with_no_covg=0; //meaning no covg on the allele we expect to have zero covg
  int num_snps_with_covg1 =0; //ignoring SNps where covg>1, likely due to repeats - this is an approx!!

  // long long total_covg_on_true_alleles=0; //2nd read in each pair
  //long long total_covg_on_error_alleles=0; //first read in each pair
  int i;
  int counts[101];
  for (i=0; i<=100; i++)
  {
    counts[i]=0;
  }
  int counts_true[1001];
  for (i=0; i<=1000; i++)
  {
    counts_true[i]=0;
  }
  
  FILE* fptr = fopen(fasta, "r");
  if (fptr==NULL)
  {
    die("Unable to open file Z%sZ. Abort.\n", fasta);
  }
  
  int num_kmers = 0;
  boolean full_entry = true;

  do
  {
    num_kmers
      = align_next_read_to_graph_and_return_node_array(fptr, max_read_length,
                                                       array_nodes, array_or,
                                                       false, &full_entry,
                                                       file_reader, seq,
                                                       kmer_window, db_graph,
                                                       -1);

    if(num_kmers > 0)
    {
      boolean too_short = false;

      long long covg_on_error_allele
        = count_reads_on_allele_in_specific_colour(array_nodes, num_kmers,
                                                   colour, &too_short);

      if(covg_on_error_allele == 0)
      {
        num_snps_with_no_covg++;
      }
      else if(covg_on_error_allele == 1)
      {
        num_snps_with_covg1++;
      }
      if(covg_on_error_allele <= 100)
      {
        counts[covg_on_error_allele] = counts[covg_on_error_allele] + 1;
      }

      num_kmers
        = align_next_read_to_graph_and_return_node_array(fptr, max_read_length,
                                                         array_nodes, array_or,
                                                         false, &full_entry,
                                                         file_reader, seq,
                                                         kmer_window, db_graph,
                                                         -1);
      
      if(num_kmers > 0)
      {
        long long covg_on_true_allele
          = count_reads_on_allele_in_specific_colour(array_nodes, num_kmers,
                                                     colour, &too_short);

        if(covg_on_true_allele <= 10000)
        {
          counts_true[covg_on_true_allele] = counts_true[covg_on_true_allele] + 1;
        }
      }
    }
  }
  while(num_kmers > 0 || !feof(fptr));

  // Finished reading, close file
  fclose(fptr);

  // Dump distribution to next line of file
  if(fout != NULL)
  {
    int j;
    fprintf(fout, "%d\t%s\t", colour, "ref_covg_distrib");

    for(j = 0; j <= 100; j++)
    {
      fprintf(fout, "%d", counts[j]);
      if(j < 100)
      {
        fprintf(fout, "\t");
      }
      else
      {
        fprintf(fout, "\n");
      }
    }

    fprintf(fout, "%d\t%s\t", colour, "alt_covg_distrib");

    for(j = 0; j <= 1000; j++)
    {
      fprintf(fout, "%d", counts_true[j]);

      if(j < 1000)
      {
        fprintf(fout, "\t");
      }
      else
      {
        fprintf(fout, "\n");
      }
    }
  }

  if(num_snps_with_no_covg + num_snps_with_covg1 == 0)
  {
    return default_seq_err;
  }
  else
  {
    printf("Num with covg 1 is %d, and with 0 is %d\n",
           num_snps_with_covg1, num_snps_with_no_covg);

    return (long double)num_snps_with_covg1 /
           (long double)(num_snps_with_no_covg + num_snps_with_covg1);
  }
}
