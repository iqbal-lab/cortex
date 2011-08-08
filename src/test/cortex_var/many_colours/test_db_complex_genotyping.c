#include <stdio.h>
#include <time.h>
#include <CUnit.h>
#include <Basic.h>
#include <test_db_complex_genotyping.h>
#include <db_complex_genotyping.h>
#include <stdlib.h>
#include <element.h>
#include <file_reader.h>
#include <simulator.h>


void test_initialise_multiplicities_of_allele_nodes_wrt_both_alleles()
{
  //create a bunch of nodes, all o which have non-zero edges in colour 0.
  //don't make any effort to create them properly - this function we are testing
  //only looks at the memory address of the node, and whether it has an edge in the right colour

  dBNode* e1 = new_element();
  e1->individual_edges[0]=1;
  dBNode* e2 = new_element();
  e2->individual_edges[0]=1;
  dBNode* e3 = new_element();
  e3->individual_edges[0]=1;
  dBNode* e4 = new_element();
  e4->individual_edges[0]=1;
  dBNode* e5 = new_element();
  e5->individual_edges[0]=1;
  dBNode* e6 = new_element();
  e6->individual_edges[0]=1;
  dBNode* e7 = new_element();
  e7->individual_edges[0]=1;
  dBNode* e8 = new_element();
  e8->individual_edges[0]=1;
  dBNode* e9 = new_element();
  e9->individual_edges[0]=1;
  dBNode* e10 = new_element();
  e10->individual_edges[0]=1;


  //first test - two branches of equal length, each with globally unique nodes
  dBNode* array1[4]={e1,e2,e3,e4};
  dBNode* array2[4]={e5,e6,e7,e8};
  VariantBranchesAndFlanks var;
  set_variant_branches_and_flanks(&var,
				  NULL, NULL, 0,
				  array1, NULL, 4,
				  array2, NULL, 4,
				  NULL, NULL, 0,
				  unknown);

  
  int mult11[4];
  int mult22[4];
  int mult12[4];
  int mult21[4];
  MultiplicitiesAndOverlapsOfBiallelicVariant mobv;
  mobv.mult11=mult11;
  mobv.mult22=mult22;
  mobv.mult12=mult12;
  mobv.mult21=mult21;
  initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, &mobv, true, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);

  int p;
  for (p=0; p<3; p++)
    {
      CU_ASSERT(mobv.mult11[p]==1);
      if (mobv.mult11[p] != 1)
	{
	  printf("Failed for p=%d. mult11[p] =%d\n", p, mobv.mult11[p]);
	}
      CU_ASSERT(mobv.mult22[p]==1);
      CU_ASSERT(mobv.mult12[p]==0);
      CU_ASSERT(mobv.mult21[p]==0);
    }


  //second test - two branches of equal length, one branch containing a repeated node
  array1[1]=e1;//so array1[0]==array1[1]
  initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, &mobv, true, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);

  CU_ASSERT(mobv.mult11[0]==2);
  CU_ASSERT(mobv.mult11[1]==2);
  CU_ASSERT(mobv.mult22[0]==1);
  CU_ASSERT(mobv.mult22[1]==1);
  CU_ASSERT(mobv.mult12[0]==0);
  CU_ASSERT(mobv.mult12[1]==0);
  CU_ASSERT(mobv.mult21[0]==0);
  CU_ASSERT(mobv.mult21[1]==0);

  for (p=2; p<4; p++)
    {
      CU_ASSERT(mobv.mult11[p]==1);
      if (mobv.mult11[p] != 1)
	{
	  printf("Failed for p=%d. mult11[p] =%d\n", p, mobv.mult11[p]);
	}
      CU_ASSERT(mobv.mult22[p]==1);
      CU_ASSERT(mobv.mult12[p]==0);
      CU_ASSERT(mobv.mult21[p]==0);
    }
  

  //third test - two branches, first shorter than second, with some nodes repeated within and between branches
  array1[1]=e2;//reset the second element in array1 to what it was before the last test
  //recall array1[4]={e1,e2,e3,e4}
  dBNode* array3[7]={e4,e1,e9,e10,e9,e9,e2};


  set_variant_branches_and_flanks(&var,
				  NULL, NULL, 0,
				  array1, NULL, 4,
				  array3, NULL, 7,
				  NULL, NULL, 0,
				  unknown);

  
  int mult3_11[4];
  int mult3_22[7];
  int mult3_12[4];
  int mult3_21[7];
  mobv.mult11=mult3_11;
  mobv.mult22=mult3_22;
  mobv.mult12=mult3_12;
  mobv.mult21=mult3_21;
  initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, &mobv, true, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);



  CU_ASSERT(mobv.mult11[0]==1);
  CU_ASSERT(mobv.mult11[1]==1);
  CU_ASSERT(mobv.mult11[2]==1);
  CU_ASSERT(mobv.mult11[3]==1);

  CU_ASSERT(mobv.mult22[0]==1);
  CU_ASSERT(mobv.mult22[1]==1);
  CU_ASSERT(mobv.mult22[2]==3);
  CU_ASSERT(mobv.mult22[3]==1);
  CU_ASSERT(mobv.mult22[4]==3);
  CU_ASSERT(mobv.mult22[5]==3);
  CU_ASSERT(mobv.mult22[6]==1);

  CU_ASSERT(mobv.mult12[0]==1);
  CU_ASSERT(mobv.mult12[1]==1);
  CU_ASSERT(mobv.mult12[2]==0);
  CU_ASSERT(mobv.mult12[3]==1);

  CU_ASSERT(mobv.mult21[0]==1);
  CU_ASSERT(mobv.mult21[1]==1);
  CU_ASSERT(mobv.mult21[2]==0);
  CU_ASSERT(mobv.mult21[3]==0);
  CU_ASSERT(mobv.mult21[4]==0);
  CU_ASSERT(mobv.mult21[5]==0);
  CU_ASSERT(mobv.mult21[6]==1);

  free(e1);
  free(e2);
  free(e3);
  free(e4);
  free(e5);
  free(e6);
  free(e7);
  free(e8);
  free(e9);
  free(e10);

}



void utility_func_test_complex_genotyping_given_two_alleles(char* fasta_allele1, char* fasta_allele2, char* fasta_genome_minus_site, char* fasta_alleles_then_genome,
							    char* fasta_allele1_then_allele2, int read_len, int kmer, int genome_size, int number_of_bits, int bucket_size,
							    int number_repeats_of_sim)
{
  if (NUMBER_OF_COLOURS<4)
    {
      printf("Need >=4 colouyrs for test_calc_log_likelihood_of_genotype_with_complex_alleles1 - recompile\n");
      exit(1);
    }

  //first set up the hash/graph
  int kmer_size = 31;
  //int number_of_bits = 15;
  //int bucket_size    = 100;
  long long bad_reads = 0; 
  long long dup_reads=0;
  int max_retries=10;
  boolean remove_duplicates_single_endedly=false; 
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  long long seq_read=0;
  long long seq_loaded=0;
  int max_read_length=2000;
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(fasta_allele1,
								     &seq_read, &seq_loaded,&bad_reads, &dup_reads, max_read_length, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,db_graph, individual_edge_array,0);

  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(fasta_allele2,
								     &seq_read, &seq_loaded,&bad_reads, &dup_reads, max_read_length, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,db_graph, individual_edge_array,0);
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(fasta_genome_minus_site, 
  								     &seq_read, &seq_loaded,&bad_reads, &dup_reads, max_read_length, 
  								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,db_graph, individual_edge_array,0);
  int max_allele_length = 90000;

  int depth;
  double seq_err_per_base;
  int colour_allele1 = 0;
  int colour_allele2 =1;
  int colour_ref_minus_site=2;
  int colour_indiv = 3;
  
  //for all of these fix read length and kmer
  read_len=50;
  kmer = 31;
  number_repeats_of_sim=100;
  

  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_allele_length,MAX_READ_NAME_LEN);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  
  
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_allele_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  
  
  //we are going to need to hold ALL the allele paths in memory at the same time
  int number_paths=3;//actiually only 2 alleles, but also one genome-minus-site will be stored in this
  dBNode*** array_of_node_arrays = (dBNode***) malloc( sizeof(dBNode**) * number_paths );
  Orientation** array_of_or_arrays = (Orientation**) malloc( sizeof(Orientation*) * number_paths );
  int* lengths_of_alleles = (int*) malloc(sizeof(int) * number_paths);
  char** array_of_allele_names = (char**) malloc( sizeof(char*) * number_paths );
  
  if ( (array_of_node_arrays==NULL) || (array_of_or_arrays==NULL) || (lengths_of_alleles==NULL) || (array_of_allele_names==NULL) )
    {
      printf("Cannot alloc arrays of arrays in test_db_complex_genotyping.c\n");
      exit(1);
    }
  
  int i;
  for (i=0; i<number_paths; i++)
    {
      array_of_node_arrays[i] = (dBNode**) malloc(sizeof(dBNode*) * max_allele_length);
      array_of_or_arrays[i]   = (Orientation*) malloc( sizeof(Orientation) * max_allele_length);
      array_of_allele_names[i]= (char*) malloc(sizeof(char) * 200 );
      
      if ( (array_of_node_arrays[i]==NULL) || (array_of_or_arrays[i]==NULL) || (array_of_allele_names==NULL) )
	{
	  printf("Cannot alloc the %d -th node and or array in print_log_liks_of_specified_set_of_genotypes_of_complex_site", i);
	  exit(1);
	}
      lengths_of_alleles[i]=0;
      array_of_allele_names[i][0]='\0';
    }



  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv = alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(max_allele_length, max_allele_length);
  if (mobv==NULL)
    {
      printf("Failed to alloc mobv\n");
      exit(1);
    }
  int* working_array_self = (int*) malloc(sizeof(int) * max_allele_length);
  int* working_array_shared = (int*) malloc(sizeof(int) * max_allele_length);
  
  if ( (mobv==NULL)||(working_array_self==NULL) || (working_array_shared==NULL))
    {
      printf("Cannot alloc all the arrays in calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site.  Give up and exit.");
      exit(1);
    }
  

  //create file reader
  int file_reader(FILE * fp, Sequence * seq, int max_allele_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      printf("new_entry must be true in hsi test function");
      exit(1);
    }
    ret =  read_sequence_from_fasta(fp,seq,max_allele_length,new_entry,full_entry,offset);
    
    return ret;
  }
  
  //end of initialisation 



  FILE* fp = fopen(fasta_alleles_then_genome, "r");
  if (fp==NULL)
    {
      printf("UNable to open %s. Exit", fasta_alleles_then_genome);
      exit(1);
    }



  int j;
  for (j=0; j<number_paths; j++)//will read allele1, then allele2, then the genome-minus-site into these arrays
    {
      int num_kmers = align_next_read_to_graph_and_return_node_array(fp, max_allele_length, 
								     array_of_node_arrays[j],
								     array_of_or_arrays[j], 
								     false, file_reader, seq, kmer_window, db_graph, -1);
    
      strcat(array_of_allele_names[j], seq->name);
      lengths_of_alleles[j]=num_kmers;
    }
  fclose(fp);
  
  VariantBranchesAndFlanks var;
  int depths[]={40,20};
  int num_depths=2;

  int p;
  GraphInfo ginfo;
  graph_info_initialise(&ginfo);
  GraphAndModelInfo model_info;

  float repeat_geometric_param_mu = 0.8;//not used in this
  // int genome_size = 554;//554 is length of one allele + rest of reference
  int num_chroms_in_expt=2;
  initialise_model_info(&model_info, &ginfo, genome_size, repeat_geometric_param_mu, 0.01, -1, num_chroms_in_expt, EachColourADiploidSample);
  boolean is_true_hom;
  for (p=0; p<num_depths; p++)
    {
      depth=depths[p];
      printf("\n********************************* New set of tests - depth = %d\n", depth);
      printf("First - tests where correct answer is het\n");
      is_true_hom=false;

      // set up the variant so it has two different alleles, one on each branch
      set_variant_branches_but_flanks_to_null(&var, 
					      array_of_node_arrays[0], 
					      array_of_or_arrays[0], 
					      lengths_of_alleles[0],
					      array_of_node_arrays[1], 
					      array_of_or_arrays[1],
					      lengths_of_alleles[1],
					      unknown);
      time_t now;
      time(&now);
      printf("%s", ctime(&now));

      printf("Branches set\n");
      initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, mobv, false, NULL, NULL);  
      time(&now);
      printf("%s", ctime(&now));

      printf("initialised multiplicities\n");
      graph_info_set_seq(model_info.ginfo, colour_indiv, genome_size*depth);
      graph_info_set_mean_readlen(model_info.ginfo, colour_indiv, read_len);

      char true_gt[]="allele1/allele2";
      // run for various sequencing error rates
      seq_err_per_base=0.01;
      model_info.seq_error_rate_per_base=seq_err_per_base;
      printf("Seq err rate is %f\n", seq_err_per_base); 
      time(&now);
      printf("%s", ctime(&now));

      printf("Call simulator\n");
      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		&var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site  fasta_allele1_then_allele2
		is_true_hom,mobv, &model_info, fasta_allele1_then_allele2, true_gt, db_graph);

      /*      
      seq_err_per_base=0.02;
      model_info.seq_error_rate_per_base=seq_err_per_base;
      printf("Seq err rate is %f\n", seq_err_per_base); 

      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		&var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		is_true_hom,mobv, &model_info, fasta_allele1_then_allele2, true_gt, db_graph);

      seq_err_per_base=0.001;
      model_info.seq_error_rate_per_base=seq_err_per_base;
      printf("Seq err rate is %f\n", seq_err_per_base); 

      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		&var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		is_true_hom,mobv, &model_info, fasta_allele1_then_allele2, true_gt, db_graph);
      
      
      // now re-run for a true hom
      printf("Now - tests where correct answer is hom1\n");
      is_true_hom=true;
      set_variant_branches_but_flanks_to_null(&var, 
					      array_of_node_arrays[0], 
					      array_of_or_arrays[0], 
					      lengths_of_alleles[0],
					      array_of_node_arrays[0], 
					      array_of_or_arrays[0],
					      lengths_of_alleles[0],
					      unknown);
      
      initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, mobv, false, NULL, NULL);
      */
      /*      
      char true_gt2[]="allele1/allele1";
      // run for various sequencing error rates
      seq_err_per_base=0.01;
      model_info.seq_error_rate_per_base=seq_err_per_base;
      printf("Seq err rate is %f\n", seq_err_per_base); 
      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		&var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		is_true_hom, mobv, &model_info, fasta_allele1_then_allele2, true_gt2, db_graph);
      
      seq_err_per_base=0.02;
      model_info.seq_error_rate_per_base=seq_err_per_base;
      printf("Seq err rate is %f\n", seq_err_per_base); 

      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		&var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		is_true_hom, mobv, &model_info, fasta_allele1_then_allele2, true_gt2, db_graph);
      seq_err_per_base=0.001;
      model_info.seq_error_rate_per_base=seq_err_per_base;
      printf("Seq err rate is %f\n", seq_err_per_base); 

      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		&var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		is_true_hom, mobv, &model_info, fasta_allele1_then_allele2, true_gt2, db_graph);
      */
    }

}


void test_calc_log_likelihood_of_genotype_with_complex_alleles1()
{

  // load a simple case. 2 alleles which do not overlap at all, or with the rest of the genome
  // load the alleles, and the rest of the genome into the hash, then run the simulator, and check the right genotype
  // happens
  utility_func_test_complex_genotyping_given_two_alleles("../data/test/pop_graph/variations/complex_genotyping/simple_allele1.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/simple_allele2.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/simple_rest_of_genome.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/simple_alleles_then_ref_minus_site.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/simple_alleles.fa", 50,31,554,15,100,100);
}




void test_calc_log_likelihood_of_genotype_with_complex_alleles2()
{

  // load a simple case. 2 alleles which do not overlap at all, or with the rest of the genome
  // load the alleles, and the rest of the genome into the hash, then run the simulator, and check the right genotype
  // happens
  utility_func_test_complex_genotyping_given_two_alleles("../data/test/pop_graph/variations/complex_genotyping/hlab_070201.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_550104.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/chr6_minus_hlab_excerpt.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_two_alleles_then_ref.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/both_hlab_alleles.fa", 50,31,89978,24,100,1);
}
