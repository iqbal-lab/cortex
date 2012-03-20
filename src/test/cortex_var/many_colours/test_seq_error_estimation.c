/*
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

#include <stdio.h>
#include <time.h>
#include <CUnit.h>
#include <Basic.h>
#include <stdlib.h>
#include <maths.h>
#include <element.h>
#include <file_reader.h>
#include <seq_error_rate_estimation.h>
#include <test_seq_error_estimation.h>

void test_estimate_seq_error_rate_for_one_colour_from_snp_allele_fasta()
{
  if (NUMBER_OF_COLOURS!=11)
    {
      printf("This test assumes 11 colours. recompile\n");
      return;
    }
  //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 10;
  int bucket_size    = 4;
  long long bad_reads = 0; 
  long long dup_reads=0;
  int max_retries=10;
  boolean remove_duplicates_single_endedly=false; 
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  int max_read_length= 200;
  long long* readlen_distrib=(long long*) malloc(sizeof(long long) * (max_read_length+1));
  long long** readlen_distrib_ptrs = (long long**) malloc(sizeof(long long*) * (max_read_length+1));
  
  if ( (readlen_distrib==NULL) || (readlen_distrib_ptrs==NULL) )
    {
      printf("Unable to malloc array to hold readlen distirbution!Exit.\n");
      exit(1);
    }
  int i;
  for (i=0; i<=max_read_length; i++)
    {
      readlen_distrib[i]=0;
      readlen_distrib_ptrs[i]=&readlen_distrib[i];
    }
  

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  long long seq_read=0;
  long long seq_loaded=0;

  //genome length 89bp, basically load 20 copies of the same genome
  load_se_and_pe_filelists_into_graph_of_specific_person_or_pop(true, false, 
								"../data/test/pop_graph/seq_error_estimation/list_sample1", 
								NULL, NULL,
								&seq_read, &seq_loaded, readlen_distrib_ptrs,
								0,false, false, false, 
								0, 33, FASTA, max_read_length, 0, db_graph);


  //initialise a graph_info object
  GraphInfo ginfo;
  graph_info_initialise(&ginfo);
  graph_info_update_mean_readlen_and_total_seq(&ginfo, 0, calculate_mean(readlen_distrib, (long long) (max_read_length+1)), seq_loaded);  
  
  long double default_err_rate = 0.01;
  estimate_seq_error_rate_from_snps_for_each_colour("../data/test/pop_graph/seq_error_estimation/colour_list_test_sample1", &ginfo, db_graph, -1, 89, default_err_rate);

  CU_ASSERT_DOUBLE_EQUAL(ginfo.seq_err[0],0.0, 0.0001);
  printf("Answer we get is %Lf and its difference from 0.5 is %Lf \n", ginfo.seq_err[0], abs(ginfo.seq_err[0]- (long double) 0.5));
}
