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

#include <CUnit.h>
#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <string.h>
#include <cmd_line.h>
#include <time.h>
#include <graph_info.h>
#include <db_differentiation.h>

void test_get_log_bayesfactor_varmodel_over_repeatmodel()
{

  if (NUMBER_OF_COLOURS<100)
    {
      printf("\nThis test requires compilation for support of >=100 colours - skipping\n");
      return;
    }

  // load a single fasta, containing a single bubble. Then run several iterations of a local function
  // that resets all the coverage values on the two branches in many colours to different scenarios.
  //Scenarios I want to test are

  // 1. All colours have both alleles with equal coverage - is repeat, so negative bayes factor
  // 2. Allele has 50% frequency, so 50% of people are hets, and 25 of each type of hom - call as a variant


  void set_coverage_on_bubble(int br1cov, int br2cov, VariantBranchesAndFlanks* var, int colour)
  {
    int i;
    for (i=0; i<var->len_one_allele; i++)
      {
	db_node_set_coverage(var->one_allele[i], individual_edge_array, colour, br1cov);
      }
    for (i=0; i<var->len_other_allele; i++)
      {
	db_node_set_coverage(var->other_allele[i], individual_edge_array, colour, br2cov);
      }
  } 


  void clean_up(dBGraph* db_graph)
  {
    void wipe_node(dBNode* node)
    {
      int i;
      for (i=0; i<100; i++)
	{
	  node->coverage[i]=0;
	}
    }
    hash_table_traverse(&wipe_node, db_graph);
  }




  //first set up the hash/graph
  int kmer_size = 13;
  int number_of_bits = 10;
  int bucket_size    = 10;
  long long bad_reads=0; 
  long long dup_reads=0;

  boolean remove_duplicates_single_endedly=false; 
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  
  
  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

  //****************************************************************************************************************
  // Then load a single bubble into the graph
  //****************************************************************************************************************
  

  

  long long seq_read=0;
  long long seq_loaded=0;
  int max_chunk_length = 30;
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/pop_graph/example1_for_testing_genotyping.allele1.fasta",&seq_read, &seq_loaded,&bad_reads, &dup_reads, max_chunk_length, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
								     db_graph, individual_edge_array, 0);
  
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/pop_graph/example1_for_testing_genotyping.allele2.fasta",&seq_read, &seq_loaded,&bad_reads, &dup_reads, max_chunk_length, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
								     db_graph, individual_edge_array, 0);


  //Now read each allele into an array of nodes, so we can use them
  

  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;

    if (new_entry!= true){
      printf("new_entry has to be true for read_next_variant_from_full_flank_file\n");
      exit(1);
    }

    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }
  
  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  int max_read_length = 30;
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
      printf("Failed to malloc kmer sliding window in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  

  int length_of_arrays=200;
  dBNode**     br1_path        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); 
  Orientation* br1_or          = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays); 
  dBNode**     br2_path        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); 
  Orientation* br2_or          = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays); 
  
  if ( (br1_path==NULL) || (br2_path==NULL) || (br1_or==NULL) || (br2_or==NULL) )
    {
      printf("Unable to alloc array for genotyping test\n");
      exit(1);
    }
  
  //end of intialisation 

  FILE* br1_fptr = fopen("../data/test/pop_graph/example1_for_testing_genotyping.allele1.fasta", "r");
  FILE* br2_fptr = fopen("../data/test/pop_graph/example1_for_testing_genotyping.allele2.fasta", "r");
  if ( (br1_fptr==NULL) || (br2_fptr==NULL) )
    {
      printf("Cannot open one of ../data/test/pop_graph/example1_for_testing_genotyping.allele1.fasta and ../data/test/pop_graph/example1_for_testing_genotyping.allele2.fasta\n");
      exit(1);
    }
  
  int br1len = align_next_read_to_graph_and_return_node_array(br1_fptr, max_read_length, br1_path, br1_or,  true, file_reader,
							       seq, kmer_window, db_graph, 0);
  int br2len = align_next_read_to_graph_and_return_node_array(br2_fptr, max_read_length, br2_path, br2_or,  true, file_reader,
							       seq, kmer_window, db_graph, 0);
  fclose(br1_fptr);
  fclose(br2_fptr);
  VariantBranchesAndFlanks var;
  var.one_allele     = br1_path;
  var.len_one_allele = br1len;
  var.other_allele     = br2_path;
  var.len_other_allele = br2len;

  int i;
  GraphInfo ginfo;
  graph_info_initialise(&ginfo);
  for (i=0; i<100; i++)
    {
      long long total_seq = 4000;
      int mean_read_len = 24; //these choices explained below
      graph_info_set_seq(&ginfo, i, total_seq);
      graph_info_set_mean_readlen(&ginfo, i, mean_read_len);
    }
  GraphAndModelInfo model_info;
  long long genome_len=100;
  double mu = 0.8; //param of geometric describing repeat copy num
  double err = 0.01;
  initialise_model_info(&model_info, &ginfo, genome_len, mu, err);
  AnnotatedPutativeVariant annovar;
  int ref_col=-1; //no reference


  clean_up(db_graph);


  //Now we are ready to start

  // 1. Repeat. Every haplotype has two copies. Lets say read length is 24, depth is 40x, k is 13. Then expected depth on a single
  //    allele is 20*(24-13+1)/24 = 10.
  //    Set up so all colours have both alleles with coverage 20

  //   Genome length? Total covg? Fix them so total covg/genome length = 40
  // say genome length=100, sequence covg = 4000;

  for (i=0; i<100; i++)
    {
      set_coverage_on_bubble(20, 20, &var, i);
    }
  
  initialise_putative_variant(&annovar, &var, BubbleCaller, &ginfo, model_info.seq_error_rate_per_base, genome_len,kmer_size,-1);
  
  double ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);

  CU_ASSERT(ret<0);//called as a repeat
  clean_up(db_graph);


  // 2.  Allele has 50% frequency, so 50% of people are hets, and 25 of each type of hom - call as a variant
  //     We tell it read_length was 12, depth was 20, so expected depth at hom site is 10.
  //     Genome length? Total covg? Fix them so total covg/genome length = 20
  //     say genome length=100, sequence covg = 2000;  

  //   Colours 0 to 49 will be het. Give both alleles covg 5 for these
  //   Colours 50 to 74 will be hom_one, covg 10 on first allele and 0 on the other
  //   Colours 75-99 will be hom_other, covg 10 on second allele and 0 on first.

  for (i=0; i<50; i++)
    {
      set_coverage_on_bubble(5,5, &var, i);
    }
  for (i=50; i<75; i++)
    {
      set_coverage_on_bubble(10,0, &var, i);
    }
  for (i=75; i<100; i++)
    {
      set_coverage_on_bubble(0,10, &var, i);
    }

  initialise_putative_variant(&annovar, &var, BubbleCaller, &ginfo, model_info.seq_error_rate_per_base, genome_len,kmer_size,ref_col);

  ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);

  CU_ASSERT(ret>0);//called as a variant

  for (i=0; i<50; i++)
    {
      CU_ASSERT(annovar.genotype[i]==het); 
    }
  for (i=50; i<75; i++)
    {
      CU_ASSERT(annovar.genotype[i]==hom_one); 
    }
  for (i=75; i<100; i++)
    {
      CU_ASSERT(annovar.genotype[i]==hom_other); 
    }
  clean_up(db_graph);







  // 3.  Allele1 has 95% frequency, so 10% of people are hets, 90% of people are hom_one, and 0.2% of people are hom_other
  //     which we treat as zero.

  //   Colours 0 to 89 will be hom_one. Give both alleles covg 5 for these
  //   Colours 90 to 99 will be het, covg 10 on first allele and 0 on the other

  for (i=0; i<90; i++)
    {
      set_coverage_on_bubble(10,0, &var, i);
    }
  for (i=90; i<100; i++)
    {
      set_coverage_on_bubble(5,5, &var, i);
    }

  initialise_putative_variant(&annovar, &var, BubbleCaller, &ginfo, model_info.seq_error_rate_per_base, genome_len,kmer_size, ref_col);
  
  ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);

  CU_ASSERT(ret>0);//called as a variant
  for (i=0; i<90; i++)
    {
      CU_ASSERT(annovar.genotype[i]==hom_one); 
    }
  for (i=90; i<100; i++)
    {
      CU_ASSERT(annovar.genotype[i]==het); 
    }
  clean_up(db_graph);





  // 4.  Allele1 has 95% frequency, so 10% of people are hets, 90% of people are hom_one, and 0.2% of people are hom_other
  //     which we treat as ONE PERSON <<<< this is the only difference with the last test


  //   Colours 0 to 89 will be hom_one. Give both alleles covg 5 for these
  //   Colours 90 to 98 will be het, covg 10 on first allele and 0 on the other
  //   Colour 99 will be hom_other, given covg 10 on second allele and 0 on first
  for (i=0; i<90; i++)
    {
      set_coverage_on_bubble(10,0, &var, i);
    }
  for (i=90; i<99; i++)
    {
      set_coverage_on_bubble(5,5, &var, i);
    }
  set_coverage_on_bubble(0,10, &var, 99);

  initialise_putative_variant(&annovar, &var, BubbleCaller, &ginfo, model_info.seq_error_rate_per_base, genome_len,kmer_size, ref_col);
  
  ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);
  CU_ASSERT(ret>0);//called as a variant
  for (i=0; i<90; i++)
    {
      CU_ASSERT(annovar.genotype[i]==hom_one); 
    }
  for (i=90; i<99; i++)
    {
      CU_ASSERT(annovar.genotype[i]==het); 
    }
  CU_ASSERT(annovar.genotype[i]==hom_other); 




  // 5.  Repeat, everybody has some covg of both alleles
  //   

  for (i=0; i<20; i++)
    {
      set_coverage_on_bubble(20,20, &var, i);
    }
  for (i=21; i<40; i++)
    {
      set_coverage_on_bubble(16,11, &var, i);
    }
  for (i=41; i<60; i++)
    {
      set_coverage_on_bubble(17,18, &var, i);
    }
  for (i=61; i<80; i++)
    {
      set_coverage_on_bubble(20,20, &var, i);
    }
  for (i=81; i<90; i++)
    {
      set_coverage_on_bubble(12,20, &var, i);
    }
  for (i=91; i<100; i++)
    {
      set_coverage_on_bubble(2,5, &var, i);
    }


  initialise_putative_variant(&annovar, &var, BubbleCaller, &ginfo, model_info.seq_error_rate_per_base, genome_len,kmer_size, ref_col);

  ret = get_log_bayesfactor_varmodel_over_repeatmodel(&annovar, &model_info);
  CU_ASSERT(ret<0);//called as a repeat





  clean_up(db_graph);


 

  //cleanup
  free_sequence(&seq);

  free(kmer_window->kmer);
  free(kmer_window);
  free(br1_path);
  free(br2_path);
  free(br1_or);
  free(br2_or);
  hash_table_free(&db_graph);

}
