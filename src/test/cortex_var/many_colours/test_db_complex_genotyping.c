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
#include <dB_graph_population.h>


void test_initialise_multiplicities_of_allele_nodes_wrt_both_alleles()
{

  if (NUMBER_OF_COLOURS<6)
    {
      printf("Must compile for >=6 colours. recompile\n");
      exit(1);
    }
  
  int working_colour1 = 2;
  int working_colour2 = 3;
  

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

  
  int mult11[4]={0,0,0,0};
  int mult22[4]={0,0,0,0};
  int mult12[4]={0,0,0,0};
  int mult21[4]={0,0,0,0};
  //  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv= alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(10,10);//no allele loinger than 10 in this test
  //reset_MultiplicitiesAndOverlapsOfBiallelicVariant(mobv);
  MultiplicitiesAndOverlapsOfBiallelicVariant mobv;
  mobv.mult11=mult11;
  mobv.mult22=mult22;
  mobv.mult12=mult12;
  mobv.mult21=mult21;
  mobv.len1=4;
  mobv.len2=4;
  improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, &mobv, true, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, 
								      working_colour1, working_colour2);

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
  reset_MultiplicitiesAndOverlapsOfBiallelicVariant(&mobv);
  improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, &mobv, true, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, working_colour1, working_colour2);

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

  
  int mult3_11[4]={0,0,0,0};
  int mult3_22[7]={0,0,0,0,0,0,0};
  int mult3_12[4]={0,0,0,0};
  int mult3_21[7]={0,0,0,0,0,0,0};
  mobv.mult11=mult3_11;
  mobv.mult22=mult3_22;
  mobv.mult12=mult3_12;
  mobv.mult21=mult3_21;
  mobv.len1 = 4;
  mobv.len2 = 7;

  improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, &mobv, true, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, working_colour1, working_colour2);



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

//you give it the fasta files for each allele (expect them to be the 1net fasta for both alleles, or the 2net fasta for both alleles)
void build_and_save_temp_binaries(char* filelist_binaries, 
				  char* first_allele_net_fasta, char* second_allele_net_fasta, char* stub,
				  int kmer, int number_of_bits, int bucket_size)
{
  long long bad_reads = 0; 
  long long dup_reads=0;
  int max_retries=10;
  boolean remove_duplicates_single_endedly=false; 
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer);

  long long seq_read=0;
  long long seq_loaded=0;
  int max_read_length=2000;
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(first_allele_net_fasta,
								     &seq_read, &seq_loaded,&bad_reads, &dup_reads, max_read_length, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,db_graph, individual_edge_array,0);
  char bin1[100];
  bin1[0]='\0';
  sprintf(bin1, "../data/tempfiles_can_be_deleted/%s_allele1_temp.ctx", stub);
  db_graph_dump_single_colour_binary_of_colour0(bin1, &db_node_check_status_not_pruned, db_graph, NULL);
  db_graph_wipe_colour(0, db_graph);


  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(second_allele_net_fasta,
								     &seq_read, &seq_loaded,&bad_reads, &dup_reads, max_read_length, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,db_graph, individual_edge_array,0);
  char bin2[100];
  bin2[0]='\0';
  sprintf(bin2, "../data/tempfiles_can_be_deleted/%s_allele2_temp.ctx", stub);
  db_graph_dump_single_colour_binary_of_colour0(bin2, &db_node_check_status_not_pruned, db_graph, NULL);
  db_graph_wipe_colour(0, db_graph);

  FILE* fp = fopen(filelist_binaries, "w");
  if (fp==NULL)
    {
      printf("Cannot open %s\n", filelist_binaries);
      exit(1);
    }
  fprintf(fp, "%s\n%s\n", bin1, bin2);
  fclose(fp);


  hash_table_free(&db_graph);
  
}


void utility_func_test_complex_genotyping_given_two_alleles(char* first_allele_name, char* second_allele_name,
							    char* fasta_allele1, char* fasta_allele2, char* fasta_genome_minus_site, char* fasta_alleles_then_genome,
							    char* fasta_allele1_then_allele2, int read_len, int kmer, int genome_size, int number_of_bits, int bucket_size,
							    int number_repeats_of_sim, boolean try_multiple_seq_errors,
							    boolean using_nets, 
							    char* first_allele_1net_fasta, char* first_allele_2net_fasta,
							    char* second_allele_1net_fasta, char* second_allele_2net_fasta )
{
  if (NUMBER_OF_COLOURS<6)
    {
      printf("Need >=6 colouyrs for test_calc_log_likelihood_of_genotype_with_complex_alleles - recompile\n");
      exit(1);
    }
  int colour_allele1 = 0;
  int colour_allele2 =1;
  int colour_ref_minus_site=2;
  int colour_indiv = 3;
  int working_colour1 = 4;
  int working_colour2 = 5;
  int working_colour_net1 = 6;
  int working_colour_net2 = 7;

  //first set up the hash/graph
  int kmer_size = kmer;
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

  //now we don't want to load the whole of the rest of the genome - just to annotate these nodes with whether they touch the rest of the genome,
  //***create local temporary new hash, load rest of genome, dump binary, then reload those nodes that are ALREADY in our allele1/2 hash

  dBGraph * temp_db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(fasta_genome_minus_site, 
  								     &seq_read, &seq_loaded,&bad_reads, &dup_reads, max_read_length, 
  								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,temp_db_graph, individual_edge_array,0);
  GraphInfo temp_db_graph_info;
  graph_info_set_seq(&temp_db_graph_info, 0, 1);//unnecessary - never used
  graph_info_set_mean_readlen(&temp_db_graph_info, 0, 1);//unnecessary - never used
  db_graph_dump_single_colour_binary_of_specified_colour("../data/tempfiles_can_be_deleted/ref_minus_genome.ctx", &db_node_condition_always_true,temp_db_graph,&temp_db_graph_info,0);
  hash_table_free(&temp_db_graph);
  int mean_r;
  long long tot_s;
  int clean_colour = 0;
  boolean only_load_kmers_already_in_graph = true;
  load_single_colour_binary_data_from_filename_into_graph("../data/tempfiles_can_be_deleted/ref_minus_genome.ctx", db_graph, &mean_r, &tot_s, false, individual_edge_array, 
							  colour_ref_minus_site, only_load_kmers_already_in_graph, clean_colour);


  int max_allele_length = 90000;

  int depth;
  double seq_err_per_base;

  //for all of these fix read length and kmer
  //read_len=50;
  //kmer = 31;
  //number_repeats_of_sim=100;
  

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



  //  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv = alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(max_allele_length, max_allele_length);
  //if (mobv==NULL)
  //  {
  //    printf("Failed to alloc mobv\n");
  //    exit(1);
  //  }
  int* working_array_self = (int*) malloc(sizeof(int) * max_allele_length);
  int* working_array_shared = (int*) malloc(sizeof(int) * max_allele_length);
  
// if ( (mobv==NULL)||(working_array_self==NULL) || (working_array_shared==NULL))
  if ( (working_array_self==NULL) || (working_array_shared==NULL))
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
  int depths[]={20};
  int num_depths=1;

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
      //printf("\n********************************* New set of tests - depth = %d\n", depth);
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
      //time_t now;
      //time(&now);
      //printf("%s", ctime(&now));

      graph_info_set_seq(model_info.ginfo, colour_indiv, genome_size*depth);
      graph_info_set_mean_readlen(model_info.ginfo, colour_indiv, read_len);

      char true_gt[100]="";
      strcat(true_gt,first_allele_name);
      strcat(true_gt,"/");
      strcat(true_gt,second_allele_name);
      // run for various sequencing error rates
      seq_err_per_base=0.01;
      model_info.seq_error_rate_per_base=seq_err_per_base;
      
      //time(&now);
      //printf("%s", ctime(&now));


      char filelist_1net_binaries[]="../data/tempfiles_can_be_deleted/filelist_hom1_1net_bins";//one binary per allele
      char filelist_2net_binaries[]="../data/tempfiles_can_be_deleted/filelist_hom1_2net_bins";



      if (using_nets==true)
	{
	  build_and_save_temp_binaries(filelist_1net_binaries, 
				       first_allele_1net_fasta, second_allele_1net_fasta, "net1", kmer, number_of_bits,  bucket_size);
	  build_and_save_temp_binaries(filelist_2net_binaries, 
				       first_allele_2net_fasta, second_allele_2net_fasta, "net2", kmer, number_of_bits,  bucket_size);
	}



      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		&var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site  fasta_allele1_then_allele2
		is_true_hom, &model_info, fasta_allele1_then_allele2, true_gt, db_graph, working_colour1, working_colour2,
		using_nets, filelist_1net_binaries, filelist_2net_binaries, working_colour_net1, working_colour_net2);

      if (try_multiple_seq_errors==true)
	{
	  seq_err_per_base=0.02;
	  model_info.seq_error_rate_per_base=seq_err_per_base;

	  
	  simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		    &var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		    is_true_hom, &model_info, fasta_allele1_then_allele2, true_gt, db_graph, working_colour1, working_colour2,
		    using_nets, filelist_1net_binaries, filelist_2net_binaries, working_colour_net1, working_colour_net2);
	  
	  seq_err_per_base=0.001;
	  model_info.seq_error_rate_per_base=seq_err_per_base;

	  
	  simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		    &var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		    is_true_hom, &model_info, fasta_allele1_then_allele2, true_gt, db_graph,working_colour1, working_colour2,
		    using_nets, filelist_1net_binaries, filelist_2net_binaries, working_colour_net1, working_colour_net2);
	}
      

      // now re-run for a true hom

      printf("Now for true hom\n");
      is_true_hom=true;
      set_variant_branches_but_flanks_to_null(&var, 
					      array_of_node_arrays[0], 
					      array_of_or_arrays[0], 
					      lengths_of_alleles[0],
					      array_of_node_arrays[0], 
					      array_of_or_arrays[0],
					      lengths_of_alleles[0],
					      unknown);


      

      char true_gt2[100]="";
      strcat(true_gt2,first_allele_name);
      strcat(true_gt2,"/");
      strcat(true_gt2,first_allele_name);

      // run for various sequencing error rates
      seq_err_per_base=0.01;
      model_info.seq_error_rate_per_base=seq_err_per_base;

      simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		&var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		is_true_hom,  &model_info, fasta_allele1_then_allele2, true_gt2, db_graph, working_colour1, working_colour2,
		using_nets, filelist_1net_binaries, filelist_1net_binaries, working_colour_net1, working_colour_net2);

      if (try_multiple_seq_errors==true)
	{      
	  seq_err_per_base=0.02;
	  model_info.seq_error_rate_per_base=seq_err_per_base;

	  
	  simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		    &var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		    is_true_hom,  &model_info, fasta_allele1_then_allele2, true_gt2, db_graph, working_colour1, working_colour2,
		    using_nets, filelist_1net_binaries, filelist_1net_binaries, working_colour_net1, working_colour_net2);
	  seq_err_per_base=0.001;
	  model_info.seq_error_rate_per_base=seq_err_per_base;

	  
	  simulator(depth, read_len, kmer, seq_err_per_base, number_repeats_of_sim, colour_indiv, colour_allele1, colour_allele2, colour_ref_minus_site,
		    &var, array_of_node_arrays[2], lengths_of_alleles[2],//these with the [2] are the genome-minus-site
		    is_true_hom,  &model_info, fasta_allele1_then_allele2, true_gt2, db_graph, working_colour1, working_colour2,
		    using_nets, filelist_1net_binaries, filelist_1net_binaries, working_colour_net1, working_colour_net2);
	}

    }

}


void test_calc_log_likelihood_of_genotype_with_complex_alleles1()
{

  // load a simple case. 2 alleles which do not overlap at all, or with the rest of the genome
  // load the alleles, and the rest of the genome into the hash, then run the simulator, and check the right genotype
  // happens
  utility_func_test_complex_genotyping_given_two_alleles("allele1", "allele2","../data/test/pop_graph/variations/complex_genotyping/simple_allele1.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/simple_allele2.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/simple_rest_of_genome.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/simple_alleles_then_ref_minus_site.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/simple_alleles.fa", 50,31,554,15,100,1, true,
							 false, NULL, NULL, NULL, NULL);
  // I have run this with 10 iterations, and it works, but is slow
}




void test_calc_log_likelihood_of_genotype_with_complex_alleles2()
{

  // load a simple case. 2 alleles which do not overlap at all, or with the rest of the genome
  // load the alleles, and the rest of the genome into the hash, then run the simulator, and check the right genotype
  // happens
  utility_func_test_complex_genotyping_given_two_alleles("hlab_070201_extended", "hlab_550104_extended", "../data/test/pop_graph/variations/complex_genotyping/hlab_070201.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_550104.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/chr6_minus_hlab_excerpt.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_two_alleles_then_ref.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/both_hlab_alleles.fa", 50,31,89978,24,100,1, false,
							 false, NULL, NULL, NULL, NULL);
}


void test_calc_log_likelihood_of_genotype_with_complex_alleles3()
{

  // load a simple case. 2 alleles which do not overlap at all, or with the rest of the genome
  // load the alleles, and the rest of the genome into the hash, then run the simulator, and check the right genotype
  // happens
  utility_func_test_complex_genotyping_given_two_alleles("hlab_070201_extended", "hlab_550104_extended", "../data/test/pop_graph/variations/complex_genotyping/hlab_070201.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_550104.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/chr6_minus_hlab_excerpt.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_two_alleles_then_ref.fa",
							 "../data/test/pop_graph/variations/complex_genotyping/both_hlab_alleles.fa", 50,31,89978,24,100,1, false,
							 true, 
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_070201_extended.single_errors.fasta", 
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_070201_extended.double_errors.fasta", 
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_550104_extended.single_errors.fasta", 
							 "../data/test/pop_graph/variations/complex_genotyping/hlab_550104_extended.double_errors.fasta");
}

/*
This test is based on my genotyping of a bad call - it's not even a bubble, 
i had a bug where the two branches do not rejoin with the right orientation.
The point is this - I want to know that when I genotype a site with the coverges below, 
with the graph_info taken from the chimp binary with which I made the calls - I want to make sure
you (=I) do not genotype lots of colours as het when they have no covg on branch1  at all!!

>var_827009_5p_flank length:1031 average_coverage:17.77 min_coverage:5 max_coverage:40 fst_coverage:15 fst_kmer:CTCTCTTGGTTGCTAATTTGTAGCAGTATGA fst_r:C fst_f:C lst_coverage:25 lst_kmer:CCTGAAGGATACAAACCCTATGATTTATCAT lst_r:A lst_f:GT 
CTCTCTTGGTTGCTAATTTGTAGCAGTATGACTCAAATGCACATAATGATAAAAATGACAAGGAGTTGGGAGTTGCAGTAGGCAAAAGTAGAACTAGGAAAAGTAATCATAGTAAAATTAAATAGGAATTCCAAATAATATAAATGTTTTTTACATTTGTGTTCAAGAAATATAAAATAAGAAAAAGGCTAAACTGTTCTTTCTAGCCTCCAAACACCATCAACCAAATTTTACTTTTATAATTTATATTTGAAAAAATGAAATTCGACCATATCTTCCCCTTAATGAGGAGTGTGTGTGTGTATGTATATCAATACACCTTGTTTTCCCCTTCTTTAGTTTCATGATTCTACCCTCAAATTAAGGTATTATTTTATTTAACAATAAATATTTATTGAGTGCCTGCCATGTTCCAGGTATTGTTCTGGTATTTGACAGACATAGGTAAGCAACACAGGAAAAGATGGCTTCCTTCTTGCACTCTTATCAGAGAAGAAAAAAAAATATCATAAAACCAATAATCATAATATCATAAAAATATTTTATGTTAAAAGGCTATTAATAATATGTAACATTTTAGAATTCTACATGGGCATAATACCATGAGGGAACTTGACCTCTTACTATATTTTATGCTTTATCTTTGAAAAAGTATACATGCCAATTGTTAAGTTCTAGGAAATCATATCACATTAATTCTCTTTGATAGCTGCAAGATTCAGATTGTTCCCATTTCAAGCCTCCTTTGGATTTCTAATTGGCTTGGCTATGCATTAACATGTAGAAAACTGTGGTCATAAATTTGTATTCAGCAAGCATTGGTAGCATGGTTCAGGCAGGTGAAAATGTATAAAGACATAAGAAAATGTGCTATGCAGGCCATTTGGAGAAATACGTCCTTGAAGATTTTTTCAAGTGGTCCCTGTGAATACAATGCATACTTTAATCTCATTTTAGTATAACTTCCTTGTTTAGATTTGTTTCCACCTGTCTGTGAATTCCTGAAGGATACAAACCCTATGATTTATCAT
>branch_827009_1 length:22 average_coverage: 1.38 min_coverage:1 max_coverage:2 fst_coverage:25 fst_kmer:CCTGAAGGATACAAACCCTATGATTTATCAT fst_r:A fst_f:GT lst_coverage:21 lst_kmer:ATTTATCATGGTATACAACATATACACTAAC lst_r:CT lst_f:C 
GGTATACAACATATACACTAAC
>branch_827009_2 length:1346 average_coverage:19.46 min_coverage:7 max_coverage:37 fst_coverage:25 fst_kmer:CCTGAAGGATACAAACCCTATGATTTATCAT fst_r:A fst_f:GT lst_coverage:21 lst_kmer:GTTAGTGTATATGTTGTATACCATGATAAAT lst_r:C lst_f:CT 
TAAACAGGCTGTCAGTTCCAAGACCAAGAGTGGCCAGATTAACTACCAATGAGAAGAAACCCTATCCTAGAGTCTGAACTTATTAATCTGTGGTCTGGGCTTTTGAGCCTTGACCCAATTGATGTCAGAAACCCTAACCAGTAAGACTTGGTCATATGCCTAGAAATGTTTCAAAGATTGCTCTTCCACACCTAGATTATTTCCAATATATTATTGGACAGCTTGTGAAACTAATCCATTTTCACTGACATGACTGGGAAATTCAGTCATAAAAAGTCTACGCAAAGATCATGTAGATATAGTTTTAAAAGGAGTAGGATAAGAAAAATGGCAGAAGGATAAGAAAAATAAGATAGGCTTCGAATGGTCTTTGTAAAGAGGTTCTATAACTATTGTTTAGTTTTTCTTGTATTAATTTCTTCCATTCATTAACAATTCAACAGTGGTTTCTATTTCTTTTTCTGCCTATTTTGTTTTGTTTTGTTTCGTCTTCTCTCTTTCTCTCTGGGATCTTCCTTCTGTTTTTATAAGAGTCTTCATTTTGTCTAAGGAATTGGCTTACAAAGATGTCCTCAATTTCTCCCATTTCATTCTCTAGTATACAGCTTGCCTTCCCAACCTTTATCACTGAAATTTGTCCGAGGTGGAAAAATAATACTTTCCTTATGGCATAAGATGAACAAGTTCTCTGACCTAAAGACTGTTGTAGTCAAGAATGTCTTTACATGCTTGGTCTTGCGAGTTCTGAGAATGCTGTCATTGAGTATGCTTGGCTGATACTCGCCCCTAACTACCATGCTTGTTGCTAACTGCTGGTACTCCACTAGGTGTACCATTGTGCCAATATGTATCTTTGTTCCAATATGTATCTTTAATTTACTAAGACTTTATTTCTGTTTAGAATGATTAAATTATAGATTTTGATTTTAGACCAATATTTTAAATTTATAAAGATAAAAAGTTTTAGGAAGCTCTTGCTGCCATTTTTCTTATCTTTCTCTTTTTAAAACTATGTCTATATGATCTTTGGCGTAGGCTTTTTATGTTTGAATTTCCCAATCATGTCAGTAAAAATGGATTAGTCTCACAAGATATCCAGTAATGTATTTGAAATAACCTACATGTGGAAGAACATTCTTTGAAATATTTCTAGGACTATGAGCAAGGCTTGCTGGTTAGGGTTTCTGATATCAGTTGGGGCAATTCTCAAAAGCCCAGACCATAGATCAATAAGCTCAGACTGTAGTATGGATCTAAAAGGCAGAGTGTATTGCCAGTTTTGGCTTCTGGTGGAAAATCCAAGCTGATAACTGAGGTTAGTGTATATGTTGTATACCATGATAAAT
>var_827009_3p_flank length:0 average_coverage:19.65 min_coverage:0 max_coverage:0 fst_coverage:21 fst_kmer:GTTAGTGTATATGTTGTATACCATGATAAAT fst_r:C fst_f:CT lst_coverage:21 lst_kmer:GTTAGTGTATATGTTGTATACCATGATAAAT lst_r:C lst_f:CT 


branch1 coverages
Covg in Colour 0:
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 1:
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 2:
6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 3:
3 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 
Covg in Colour 4:
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 5:
4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 6:
2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 7:
3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 8:
2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 9:
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in Colour 10:
4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
branch2 coverages
Covg in Colour 0:
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Covg in Colour 1:
0 0 0 0 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 4 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 5 5 5 4 4 4 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 5 5 5 5 5 5 5 5 5 3 3 3 3 3 3 3 4 4 4 4 4 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 2 3 3 3 4 4 4 4 5 5 5 5 5 4 4 4 5 5 5 5 5 4 3 3 3 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 3 3 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 2 2 3 3 3 3 3 3 3 3 2 2 2 3 3 3 2 2 3 3 3 3 3 2 2 2 2 2 2 2 2 3 3 3 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 0 0 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 4 4 4 4 5 5 5 4 4 4 3 3 3 3 3 3 3 3 2 2 2 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 3 3 4 5 6 6 6 6 8 8 8 8 7 7 6 6 6 6 6 6 5 4 4 3 3 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 3 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 2 2 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 2 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 2 2 2 2 3 3 3 3 3 3 3 4 6 7 7 7 6 7 7 7 7 6 6 6 6 5 5 5 5 5 5 5 4 2 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 3 3 3 3 3 4 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Covg in Colour 2:
6 6 6 6 4 3 3 3 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 3 3 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 1 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 2 2 3 3 3 3 3 3 3 4 3 3 3 3 3 3 3 3 3 3 3 2 3 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 1 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 4 2 2 2 1 1 1 1 1 1 1 1 3 3 3 4 4 4 4 4 3 3 3 4 5 5 5 5 5 5 5 5 5 3 3 3 2 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 4 4 4 4 3 3 3 2 2 2 2 2 2 2 2 2 2 1 1 2 2 1 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2 3 5 5 6 6 6 6 6 6 5 5 5 5 5 5 4 4 5 5 5 5 4 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 3 3 3 3 3 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 3 3 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 3 3 4 3 3 3 3 3 3 2 2 2 2 3 3 4 4 4 5 5 5 4 4 3 3 3 3 3 3 3 4 5 5 4 4 4 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 5 5 4 4 4 3 3 3 3 3 3 3 3 2 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 5 5 5 5 5 4 4 3 3 3 3 3 3 3 4 4 4 4 4 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 5 5 4 4 4 3 3 3 3 4 3 3 3 3 4 3 3 3 3 3 3 3 3 3 3 4 3 3 3 3 2 2 2 2 2 1 1 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 2 2 3 3 3 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 4 4 4 3 4 4 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 3 2 3 3 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 2 2 4 4 4 4 5 5 5 5 4 4 4 4 4 4 4 4 4 5 5 4 4 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 3 5 4 4 4 5 5 5 6 6 6 6 6 7 8 8 10 10 10 10 9 8 6 6 6 6 5 5 5 4 4 4 4 4 3 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 
Covg in Colour 3:
3 2 2 2 2 2 1 2 2 2 2 2 1 1 1 1 0 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 0 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 1 2 3 3 3 4 4 4 4 4 4 4 4 4 6 6 5 6 6 6 6 6 5 4 4 4 3 3 3 3 3 3 3 3 3 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 3 3 1 1 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 1 1 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 2 2 2 2 2 3 3 3 2 2 2 2 2 3 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 3 4 3 1 2 2 2 2 2 2 3 3 3 4 3 3 3 3 3 3 3 3 2 2 2 1 1 1 1 2 2 1 1 1 1 1 2 2 2 2 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 4 4 3 3 3 3 2 2 2 2 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 3 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 3 3 3 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 2 2 3 3 3 3 4 4 4 3 3 3 4 4 4 5 4 5 5 5 5 5 5 4 3 3 3 2 2 2 2 2 2 2 2 2 1 1 0 0 0 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 4 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 2 2 3 2 3 3 4 4 4 4 5 5 5 5 5 6 6 6 6 6 6 6 6 5 5 4 5 4 5 5 5 4 4 4 4 4 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 4 4 4 4 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 3 3 3 
Covg in Colour 4:
1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 2 2 2 2 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 2 2 2 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 6 6 6 5 5 4 3 2 2 2 2 2 2 2 2 3 3 3 3 3 3 2 2 2 2 2 2 3 4 4 4 4 4 4 4 5 4 4 4 4 4 4 3 4 4 4 4 4 3 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 5 5 4 4 5 5 5 5 5 5 5 5 5 3 3 3 3 4 4 4 4 2 2 3 3 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 2 2 2 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 4 3 3 3 2 2 2 2 1 1 2 2 3 3 3 3 3 3 3 3 3 2 2 3 3 4 4 4 4 4 4 3 3 3 3 3 3 3 3 5 5 5 5 5 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 3 3 3 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 5 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 2 2 2 2 2 2 2 3 3 1 1 1 1 2 3 3 2 2 2 2 2 2 2 2 2 2 3 4 4 4 4 4 4 5 5 4 4 4 4 5 5 6 6 6 6 6 6 5 4 4 4 4 4 4 3 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 1 1 1 1 1 1 1 1 2 2 2 3 3 4 4 4 4 4 4 3 3 3 3 3 4 4 4 4 4 3 4 4 3 3 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 2 2 2 2 2 3 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 4 4 4 4 4 3 3 3 3 3 3 3 3 4 3 3 3 3 3 3 3 1 1 1 1 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 5 5 5 5 4 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 2 4 4 4 3 3 3 3 3 3 3 3 3 3 2 4 5 5 5 5 5 5 3 3 3 3 3 3 3 3 4 4 5 7 7 7 5 4 4 4 4 4 4 5 5 5 5 5 5 4 4 4 4 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 4 5 5 6 6 6 6 6 5 5 5 5 5 5 5 5 5 4 4 6 6 4 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 1 1 2 2 2 2 2 3 3 3 3 4 4 4 4 4 4 4 4 4 4 3 3 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 1 1 1 1 1 2 2 2 2 2 2 
Covg in Colour 5:
4 4 4 4 4 4 4 4 4 4 4 4 2 3 4 3 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 2 2 2 3 3 4 4 4 5 5 5 5 5 6 6 6 7 8 8 8 7 6 6 6 5 5 4 4 4 3 4 4 4 4 3 3 3 2 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 2 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 4 3 3 4 4 4 4 4 4 4 3 2 2 2 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 2 2 3 3 3 3 2 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 3 3 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 3 2 2 2 2 3 4 4 4 5 5 5 5 5 5 5 5 5 3 3 6 6 6 5 5 6 6 5 6 6 5 5 5 5 5 5 5 5 5 5 5 3 3 3 3 3 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 4 4 5 5 5 5 6 6 6 6 6 6 6 6 5 5 5 5 5 5 4 3 3 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 1 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 3 3 3 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 2 3 3 4 5 5 5 5 4 4 4 4 4 5 5 5 5 5 5 5 5 4 3 3 2 1 1 1 2 3 3 3 2 2 3 3 3 3 3 3 6 6 6 6 6 6 7 7 7 6 5 5 5 6 6 5 5 5 5 5 5 2 2 2 2 2 2 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 2 2 2 3 4 4 4 4 4 4 4 4 5 4 4 4 4 4 4 4 3 3 3 3 4 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 1 1 1 1 1 1 1 0 0 0 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 
Covg in Colour 6:
2 2 2 2 2 2 1 1 1 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 2 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 2 2 2 2 3 3 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 2 2 2 2 2 2 3 3 3 3 4 4 4 4 4 4 4 4 4 3 3 2 2 2 2 2 2 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 2 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 3 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 4 4 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 3 3 3 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 6 6 6 5 4 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 1 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 5 5 5 6 6 6 6 6 7 7 6 6 6 6 6 6 6 6 6 6 6 4 4 4 3 3 3 2 2 1 1 1 2 2 2 2 2 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 2 2 2 2 2 2 3 3 4 4 4 4 4 5 5 5 5 5 4 4 4 4 4 4 4 4 4 3 4 3 3 3 3 3 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 4 5 5 5 5 4 4 3 3 3 3 3 3 3 3 3 3 3 1 1 1 0 1 1 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 4 4 4 4 3 3 2 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 2 3 3 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 1 2 2 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 2 2 2 2 
Covg in Colour 7:
3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 3 3 3 4 4 4 4 4 4 4 4 4 4 4 6 6 6 5 4 5 5 4 4 5 4 4 4 3 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 1 1 2 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 3 3 3 3 3 2 1 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 3 3 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 2 2 2 2 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 3 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 3 3 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 3 2 2 2 2 1 1 1 1 1 2 3 3 3 3 3 4 4 4 4 4 3 4 4 4 4 4 4 5 5 5 4 3 4 4 4 5 4 4 4 4 4 4 3 3 3 3 3 4 3 3 3 4 4 3 3 3 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 1 1 1 1 1 2 4 4 4 4 4 4 4 4 5 4 5 5 6 6 6 6 6 6 6 6 6 4 4 5 5 5 5 5 5 4 5 4 4 3 3 3 3 4 4 4 4 3 3 3 2 2 3 3 3 3 3 2 2 2 3 3 3 4 3 3 3 3 3 3 3 3 3 2 2 1 1 2 2 3 4 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 2 2 2 1 2 2 2 2 2 4 4 4 3 3 3 3 4 4 4 4 4 4 6 6 6 5 4 4 4 4 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 2 2 2 2 2 3 4 5 5 5 4 4 4 5 5 6 6 6 6 6 6 6 6 6 6 6 5 4 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 2 2 2 2 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 3 3 3 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 0 1 3 3 3 3 3 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 5 3 3 3 3 3 3 3 3 3 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 4 5 5 5 5 4 4 4 4 4 4 5 5 5 5 5 5 3 3 3 3 2 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 4 4 4 4 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 2 2 2 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 3 3 2 2 2 2 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 3 3 3 3 3 3 3 4 3 3 4 4 4 4 4 4 4 4 4 4 3 2 2 3 3 3 3 3 2 2 2 1 1 1 1 1 2 2 2 2 2 2 4 4 3 5 5 5 5 5 6 6 7 7 7 7 9 8 8 8 8 8 8 6 6 6 5 5 5 6 6 5 5 4 4 4 4 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 6 6 7 7 7 7 7 7 7 7 6 8 7 8 8 7 7 6 6 6 6 3 3 3 3 3 3 3 3 3 3 3 1 1 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 1 1 1 1 2 2 2 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 3 3 3 3 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 2 2 2 2 2 2 3 3 3 4 4 4 4 4 4 4 3 3 3 3 3 2 2 2 2 2 2 1 1 1 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 1 1 1 2 2 2 2 2 2 3 2 2 2 2 2 1 1 1 1 1 1 2 2 2 1 1 
Covg in Colour 8:
2 4 4 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 4 4 4 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 4 4 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 3 3 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 2 2 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 4 4 4 3 3 3 3 3 3 3 6 6 7 7 7 7 7 7 8 8 8 4 4 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 3 3 3 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 3 3 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 5 5 5 5 6 6 6 4 4 4 5 5 6 6 5 5 5 5 7 8 8 7 7 7 7 7 7 7 6 7 7 6 6 6 6 6 8 8 8 6 5 5 5 5 5 5 3 3 3 3 3 3 3 5 5 5 5 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 4 4 3 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 3 3 3 3 3 3 2 2 2 2 2 2 2 3 4 4 4 4 4 4 4 2 2 2 2 2 3 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 3 3 3 3 3 5 5 8 8 8 8 8 8 8 9 9 9 8 7 7 7 7 7 7 7 7 6 6 4 4 4 5 5 4 4 1 1 1 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 0 0 0 1 1 1 1 2 2 2 2 2 4 4 4 4 5 5 5 5 5 5 5 5 4 4 4 4 3 3 3 3 3 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 4 4 4 4 4 3 3 3 3 3 3 3 4 4 4 4 4 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 3 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 3 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 0 0 0 0 1 1 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 5 6 2 2 2 2 2 2 3 3 3 4 4 4 4 4 4 4 3 3 3 3 2 2 3 3 3 3 3 2 2 3 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 5 5 3 3 6 6 6 5 5 5 5 5 5 5 5 4 4 4 4 4 4 3 3 3 3 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 2 2 2 3 3 3 3 3 5 5 
Covg in Colour 9:
0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 3 3 3 2 2 2 2 2 3 3 3 3 3 3 5 5 7 6 6 6 6 5 5 5 5 5 5 4 5 5 5 5 5 5 5 3 3 2 2 2 2 2 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 2 3 3 3 3 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 3 4 4 5 5 5 5 5 5 5 5 3 3 3 3 3 3 3 3 3 3 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 0 2 2 2 3 3 3 3 3 3 2 2 2 2 4 4 4 5 5 5 5 5 3 3 4 4 4 4 4 4 4 4 4 4 4 2 2 2 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 2 2 2 2 2 2 2 1 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 2 2 2 2 3 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 5 4 4 4 4 3 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 4 4 4 5 4 5 4 4 4 4 4 3 3 3 3 3 2 2 2 2 2 2 2 2 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 3 4 4 4 4 3 4 4 4 5 4 4 4 4 4 4 7 7 8 8 8 7 6 6 6 6 6 5 5 5 4 4 4 4 4 4 4 1 1 0 0 0 0 0 0 0 0 0 1 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 1 2 2 2 3 3 3 2 2 2 2 2 3 3 3 3 3 3 2 2 2 3 3 3 3 2 2 2 3 3 4 4 4 3 3 4 4 4 4 4 4 4 3 3 3 3 3 3 3 2 2 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 2 2 1 1 1 2 3 4 3 3 2 3 3 3 3 3 3 3 3 3 4 5 5 5 5 5 4 4 3 3 3 3 2 2 2 3 3 3 3 3 4 3 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 2 2 2 2 1 1 1 2 2 2 1 1 1 1 2 3 3 3 3 4 5 5 5 5 4 4 4 4 5 5 5 5 5 5 5 4 3 3 3 3 2 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 2 2 2 2 2 2 2 2 3 4 5 5 5 5 5 5 4 4 4 4 5 3 3 4 4 4 4 4 4 3 2 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 2 3 3 3 4 4 4 4 5 5 4 4 4 4 4 5 5 5 5 5 5 4 3 2 3 2 2 2 3 5 5 5 5 5 5 6 5 5 5 5 5 6 6 6 6 5 5 6 6 5 3 3 3 3 3 3 2 2 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 2 2 3 3 4 4 3 4 4 4 4 4 4 4 4 4 5 5 5 5 5 4 4 3 3 2 2 2 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 3 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 4 4 4 4 4 4 3 3 4 4 4 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 2 2 2 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 3 2 4 4 3 3 3 3 
Covg in Colour 10:
4 4 4 3 3 3 4 4 4 3 3 3 3 4 4 4 4 4 3 3 3 3 3 3 3 3 3 2 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 4 4 4 4 4 4 3 3 3 3 3 3 3 1 2 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 4 4 4 4 3 3 3 3 4 4 3 3 3 3 3 3 3 2 2 2 2 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 3 3 3 3 3 3 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 2 3 3 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 4 5 5 5 4 4 4 3 2 3 3 3 3 3 3 3 3 3 3 3 3 2 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 1 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 6 5 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 4 4 4 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 2 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 3 3 3 2 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 4 3 3 3 4 4 4 4 5 5 6 5 5 5 5 6 7 7 6 6 6 5 5 5 5 4 4 4 5 4 4 3 3 3 3 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 4 4 5 6 6 6 6 6 7 7 6 6 6 6 5 5 5 5 5 5 5 4 4 3 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 3 3 4 4 4 3 3 3 3 3 3 3 2 2 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 4 4 4 4 4 5 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 4 4 4 3 3 3 3 3 3 3 3 3 4 4 4 3 4 4 5 6 5 5 5 5 5 6 6 7 7 7 7 7 7 6 6 6 6 5 5 4 3 3 2 2 2 2 1 1 0 0 0 0 0 0 1 1 1 2 2 2 2 2 2 2 3 4 4 4 4 4 4 4 4 4 4 3 3 3 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 1 2 2 2 2 3 3 3 3 3 4 4 4 3 3 3 3 3 4 4 4 3 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 2 2 2 2 2 2 2 


 */

void regression_test_1_single_bubble_call_one_allele_shorter_than_k_one_very_long()
{
  if (NUMBER_OF_COLOURS<11)
    {
      printf("This test requires >=11 colours, skipping\n");
      return;
    }

  VariantBranchesAndFlanks var;
  dBNode* branch1[22];
  dBNode* branch2[134];
  int i,j;
  for (i=0; i<22; i++)
    {
      branch1[i]=new_element();
    }
  for (i=0; i<134; i++)
    {
      branch2[i]=new_element();
    }
  //int widths[] = { [0 ... 9] = 1, [10 ... 99] = 2, [100] = 3 };
  int br1_covg_0[]={[0]=1, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[0]=br1_covg_0[j];
    }
  int br1_covg_1[]={[0 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[1]=br1_covg_1[j];
    }
  int br1_covg_2[]={[0]=6, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[2]=br1_covg_2[j];
    }

  int br1_covg_3[]={3,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,1};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[3]=br1_covg_3[j];
    }

  int br1_covg_4[]={[0]=1, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[4]=br1_covg_4[j];
    }


  int br1_covg_5[]={[0]=4, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[5]=br1_covg_5[j];
    }

  int br1_covg_6[]={[0]=2, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[6]=br1_covg_6[j];
    }

  int br1_covg_7[]={[0]=3, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[7]=br1_covg_7[j];
    }

  int br1_covg_8[]={[0]=2, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[8]=br1_covg_8[j];
    }
  //colour 9 all zeroes, no need to do anything
  
  int br1_covg_10[]={[0]=4, [1 ... 21]=0};
  for (j=0; j<22; j++)
    {
      branch1[j]->coverage[10]=br1_covg_10[j];
    }


  //branch2 coverages
  int br2_covg_0[]={[0 ... 133]=1};
  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[0]=br2_covg_0[j];
    }

  int br2_covg_1[]={0,0,0,0,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,2,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4};
  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[1]=br2_covg_1[j];
    }


  int br2_covg_2[]={6,6,6,6,4,3,3,3,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1};
  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[2]=br2_covg_2[j];
    }


  int br2_covg_3[]={3,2,2,2,2,2,1,2,2,2,2,2,1,1,1,1,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,0,0,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[3]=br2_covg_3[j];
    }



  int br2_covg_4[]={1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,2,2,2,3,3,3,3,3,3,4,4,4,4,3,3,3,3,3,3,3,2,2,2,2,1,1,1,1,1,2,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[4]=br2_covg_4[j];
    }


  int br2_covg_5[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,2,2,2,2,3,4,4,4,5,5,5,5,5,5,5,5,5,3,3,6,6,6,5,5,6,6,5,6,6,5,5,5,5,5,5,5,5,5,5,5,3,3,3,3,3,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[5]=br2_covg_5[j];
    }


  int br2_covg_6[]={2,2,2,2,2,2,1,1,1,2,2,2,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[6]=br2_covg_6[j];
    }



  int br2_covg_7[]={3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,6,6,6,5,4,5,5,4,4,5,4,4,4,3,3,3,3,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,1,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,3,3,3,3,3,2,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[7]=br2_covg_7[j];
    }


  int br2_covg_8[]={2,4,4,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,0,0,0,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,4,4,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,4,4,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2};


  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[8]=br2_covg_8[j];
    }


  int br2_covg_9[]={0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,2,2,2,2,2,1,1,1,2,2,2,2,2,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2};

  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[9]=br2_covg_9[j];
    }



  int br2_covg_10[]={4,4,4,3,3,3,4,4,4,3,3,3,3,4,4,4,4,4,3,3,3,3,3,3,3,3,3,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,4,4,4,4,4,4,3,3,3,3,3,3,3,1,2,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,1,1,1,1,1,2,2,2,3,2,2,2,2,2,2,4,4};


  for (j=0; j<134; j++)
    {
      branch2[j]->coverage[10]=br2_covg_10[j];
    }


  /*
    output from making the chimp calls:

SUMMARY:
Colour:	MeanReadLen	TotalSeq
0	0	2912760135
1	50	18286122352
2	50	16816361244
3	50	18039181209
4	50	15879192506
5	50	17729089947
6	50	15750659112
7	50	26196361173
8	52	20202087523
9	50	18907785783
10	50	16870486574
****************************************

   */

  GraphInfo ginfo;
  graph_info_initialise(&ginfo);
  graph_info_set_seq(&ginfo, 0, 2912760135);
  graph_info_set_seq(&ginfo, 1, 18286122352 );
  graph_info_set_seq(&ginfo, 2, 16816361244);
  graph_info_set_seq(&ginfo, 3, 18039181209);
  graph_info_set_seq(&ginfo, 4, 15879192506);
  graph_info_set_seq(&ginfo, 5, 17729089947);
  graph_info_set_seq(&ginfo, 6, 15750659112);
  graph_info_set_seq(&ginfo, 7, 26196361173);
  graph_info_set_seq(&ginfo, 8, 20202087523);
  graph_info_set_seq(&ginfo, 9, 18907785783);
  graph_info_set_seq(&ginfo, 10, 16870486574);
  graph_info_set_mean_readlen(&ginfo, 0, 0);
  graph_info_set_mean_readlen(&ginfo, 1, 50);
  graph_info_set_mean_readlen(&ginfo, 2, 50);
  graph_info_set_mean_readlen(&ginfo, 3, 50);
  graph_info_set_mean_readlen(&ginfo, 4, 50);
  graph_info_set_mean_readlen(&ginfo, 5, 50);
  graph_info_set_mean_readlen(&ginfo, 6, 50);
  graph_info_set_mean_readlen(&ginfo, 7, 50);
  graph_info_set_mean_readlen(&ginfo, 8, 52);
  graph_info_set_mean_readlen(&ginfo, 9, 50);
  graph_info_set_mean_readlen(&ginfo, 10, 50);
 

  GraphAndModelInfo model_info;
  double mu=0.8;
  double seq_err_rate_per_base=0.01;
  int ref_colour=0;
  int num_chroms=20; 
  long long genome_len = 3000000000;
  initialise_model_info(&model_info, &ginfo, genome_len, mu, seq_err_rate_per_base, 
			ref_colour, num_chroms, EachColourADiploidSampleExceptTheRefColour);

  var.one_allele       = branch1;
  var.len_one_allele   = 22;
  var.other_allele     = branch2;
  var.len_other_allele = 134;

  AnnotatedPutativeVariant annovar;
  initialise_putative_variant(&annovar, &model_info, &var, BubbleCaller, 31, AssumeUncleaned, NULL, NULL, NULL );//last 3 arguments (gwp, db_graph and little db graph) are not used except for PD genotying.


  // Since none of the colours except colour 3 has any coverage AT ALL on branch1, I simply
  // cannot accept a genotype call which is het or hom_one for those colours

  CU_ASSERT(annovar.genotype[0]!=hom_one);
  CU_ASSERT(annovar.genotype[1]!=hom_one);
  CU_ASSERT(annovar.genotype[2]!=hom_one);
  //leaving out colour 3
  CU_ASSERT(annovar.genotype[4]!=hom_one);
  CU_ASSERT(annovar.genotype[5]!=hom_one);
  CU_ASSERT(annovar.genotype[6]!=hom_one);
  CU_ASSERT(annovar.genotype[7]!=hom_one);
  CU_ASSERT(annovar.genotype[8]!=hom_one);
  CU_ASSERT(annovar.genotype[9]!=hom_one);
  CU_ASSERT(annovar.genotype[10]!=hom_one);
  
  //hom one is least likely
  CU_ASSERT(annovar.gen_log_lh[10].log_lh[0]<annovar.gen_log_lh[10].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[10].log_lh[0]<annovar.gen_log_lh[10].log_lh[2]);

  CU_ASSERT(annovar.gen_log_lh[1].log_lh[0]<annovar.gen_log_lh[1].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[1].log_lh[0]<annovar.gen_log_lh[1].log_lh[2]);

  CU_ASSERT(annovar.gen_log_lh[2].log_lh[0]<annovar.gen_log_lh[2].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[2].log_lh[0]<annovar.gen_log_lh[2].log_lh[2]);

  CU_ASSERT(annovar.gen_log_lh[3].log_lh[0]<annovar.gen_log_lh[3].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[3].log_lh[0]<annovar.gen_log_lh[3].log_lh[2]);

  CU_ASSERT(annovar.gen_log_lh[4].log_lh[0]<annovar.gen_log_lh[4].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[4].log_lh[0]<annovar.gen_log_lh[4].log_lh[2]);

  CU_ASSERT(annovar.gen_log_lh[5].log_lh[0]<annovar.gen_log_lh[5].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[5].log_lh[0]<annovar.gen_log_lh[5].log_lh[2]);

  CU_ASSERT(annovar.gen_log_lh[6].log_lh[0]<annovar.gen_log_lh[6].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[6].log_lh[0]<annovar.gen_log_lh[6].log_lh[2]);

  CU_ASSERT(annovar.gen_log_lh[7].log_lh[0]<annovar.gen_log_lh[7].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[7].log_lh[0]<annovar.gen_log_lh[7].log_lh[2]);

  CU_ASSERT(annovar.gen_log_lh[8].log_lh[0]<annovar.gen_log_lh[8].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[8].log_lh[0]<annovar.gen_log_lh[8].log_lh[2]);

  CU_ASSERT(annovar.gen_log_lh[9].log_lh[0]<annovar.gen_log_lh[9].log_lh[1]);
  CU_ASSERT(annovar.gen_log_lh[9].log_lh[0]<annovar.gen_log_lh[9].log_lh[2]);




  //cleanup

  for (i=0; i<22; i++)
    {
      free(branch1[i]);
    }
  for (i=0; i<134; i++)
    {
      free(branch2[i]);
    }

}



void regression_test_2_genotyping_of_PD_SNP_call()
{
  int kmer_size = 55;

  if (NUMBER_OF_COLOURS<2)
    {
      printf("This test requires >=2 colours, skipping\n");
      return;
    }

  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      //offset = kmer_size;
      printf("new_entry must be true in hsi test function");
      exit(1);
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }



  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------

  int max_read_length=400;

  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,max_read_length);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. Exit.\n");
      exit(1);
    }
  

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  
  
  //end of intialisation 


  /*

One of the variants called on chr22 on NA12878

>var_18_5p_flank length:165 average_coverage: 5.00 min_coverage:3 max_coverage:11 mode_coverage: 4 percent_nodes_with_modal_covg: 36.00 percent_novel:  0.00 fst_coverage:11 fst_kmer:AATACTATTTTTAGGCCAGGCATGGTGGCTCATACCTGTAATCCCAACATTTTGG fst_r:A fst_f:G lst_coverage:6 lst_kmer:TAGTGAGACCTTTCCTTTATTAAATAAATAAATAGGATGGGCACTGTGGCTCATA lst_r:T lst_f:T 
AATACTATTTTTAGGCCAGGCATGGTGGCTCATACCTGTAATCCCAACATTTTGGGAGGCCAAGTTTGGAGACTCATTTGAGTCCAGGAGTTGACCAGCCTGGGCAACATAGTGAGACCTTTCCTTTATTAAATAAATAAATAGGATGGGCACTGTGGCTCATAT
>var_18_trusted_branch length:56 average_coverage: 1.00 min_coverage:0 max_coverage:19 mode_coverage: 0 percent_nodes_with_modal_covg: 92.00 percent_novel:  0.00 fst_coverage:0 fst_kmer:GTGAGACCTTTCCTTTATTAAATAAATAAATAGGATGGGCACTGTGGCTCATATC fst_r: fst_f: lst_coverage:8 lst_kmer:TGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT lst_r:C lst_f:T covgs of trusted not variant nodes:  0 0 0 0 0 0 0 0 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 19 0 0 0 0 0 0 0 0 16 0 0 0 0 0 number_of_such_nodes: 55
CTGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT
>var_18_variant_branch length:56 average_coverage: 6.00 min_coverage:4 max_coverage:9 mode_coverage: 7 percent_nodes_with_modal_covg: 46.00 percent_novel: 98.00 fst_coverage:6 fst_kmer:AGTGAGACCTTTCCTTTATTAAATAAATAAATAGGATGGGCACTGTGGCTCATAT fst_r:A fst_f:G lst_coverage:8 lst_kmer:TGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT lst_r:C lst_f:T covgs of variant not trusted nodes:  7 7 5 8 5 8 8 7 8 8 7 7 7 5 7 4 7 6 8 7 9 7 9 7 7 7 7 7 6 7 8 7 6 9 7 7 7 7 8 5 8 7 4 7 8 6 7 7 7 8 6 8 6 5 8 number_of_such_nodes: 55
GTGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT
>var_18_3p_flank length:45 average_coverage: 4.00 min_coverage:1 max_coverage:8 mode_coverage: 2 percent_nodes_with_modal_covg: 35.00 percent_novel:  0.00 fst_coverage:8 fst_kmer:TGTAATCCCAGCATTTTGGGTTGCCAAGGCAGGAGGCTTGCTTGAGCCCAGGAGT fst_r:C fst_f:T lst_coverage:1 lst_kmer:GCCCAGGAGTTTGAGACCAGCCTGGACAGCATAGCAAGACTCCATCTCTACAAAT lst_r:T lst_f: 
TTGAGACCAGCCTGGACAGCATAGCAAGACTCCATCTCTACAAAT
var_18 - extra information

branch1 coverages
Mult in  hum ref:
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
Covg in indiv:
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16 16 19 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 
branch2 coverages
Mult in  hum ref:
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
Covg in indiv:
6 6 6 7 7 7 7 7 7 7 8 8 9 9 9 8 8 8 8 8 8 7 5 6 6 6 5 5 4 4 5 5 6 7 8 8 8 7 7 7 7 7 7 7 7 7 7 7 7 8 7 7 7 7 7 8 
  */




  //first set up the hash/graph

  int number_of_bits = 18;
  int bucket_size    = 100;
  long long bad_reads = 0; 
  long long dup_reads=0;
  int max_retries=10;
  boolean remove_duplicates_single_endedly=false; 
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);


  long long seq_read=0;
  long long seq_loaded=0;
  int max_read_len = 300;

  //I load this into the graph so the kmers are there, but then I am going to just create a car object with the covgs I want

  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/pop_graph/variations/complex_genotyping/pd_example1.fasta",
								     &seq_read, &seq_loaded,
								     &bad_reads, &dup_reads, 
								     max_read_len, 
								     remove_duplicates_single_endedly, 
								     break_homopolymers, 
								     homopolymer_cutoff,db_graph,individual_edge_array,0);

  //for second var/test
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/pop_graph/variations/complex_genotyping/pd_example2.fasta",
								     &seq_read, &seq_loaded,
								     &bad_reads, &dup_reads, 
								     max_read_len, 
								     remove_duplicates_single_endedly, 
								     break_homopolymers, 
								     homopolymer_cutoff,db_graph,individual_edge_array,0);



  VariantBranchesAndFlanks var;
  dBNode** branch1 = (dBNode**) malloc(sizeof(dBNode*) * max_read_length);
  dBNode** branch2 = (dBNode**) malloc(sizeof(dBNode*) * max_read_length);
  Orientation* branch1_o = (Orientation*) malloc(sizeof(Orientation)*max_read_length);
  Orientation* branch2_o = (Orientation*) malloc(sizeof(Orientation)*max_read_length);

  if (  (branch1==NULL) || (branch2==NULL) || (branch1_o==NULL) || (branch2_o==NULL) )
    {
      printf("Unable to malloc branch arrays in test\n");
      exit(1);
    }


  FILE* fp = fopen("../data/test/pop_graph/variations/complex_genotyping/pd_example1_just_alleles.fasta", "r");
  if (fp==NULL)
    {
      printf("Unable to open test file ../data/test/pop_graph/variations/complex_genotyping/pd_example1_just_alleles.fasta\n");
      exit(1);
    }
  int num_kmers_br1 = align_next_read_to_graph_and_return_node_array(fp, max_read_length, branch1, branch1_o, false, &file_reader_fasta, 
								      seq, kmer_window, db_graph, 0);

  int num_kmers_br2 = align_next_read_to_graph_and_return_node_array(fp, max_read_length, branch2, branch2_o, false, &file_reader_fasta, 
								      seq, kmer_window, db_graph, 0);


  fclose(fp);
  int i,j;

  int br1_covg_0[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  for (j=0; j<56; j++)
    {
      branch1[j]->coverage[0]=br1_covg_0[j];
    }
  int br1_covg_1[]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,16,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8};

  for (j=0; j<56; j++)
    {
      branch1[j]->coverage[1]=br1_covg_1[j];
    }



  int br2_covg_0[]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  for (j=0; j<56; j++)
    {
      branch2[j]->coverage[0]=br2_covg_0[j];
    }
  int br2_covg_1[]= {6,6,6,7,7,7,7,7,7,7,8,8,9,9,9,8,8,8,8,8,8,7,5,6,6,6,5,5,4,4,5,5,6,7,8,8,8,7,7,7,7,7,7,7,7,7,7,7,7,8,7,7,7,7,7,8 };

  for (j=0; j<56; j++)
    {
      branch2[j]->coverage[1]=br2_covg_1[j];
    }
  
  /*

The above site was called by the PD caller on NA12878 with k=55 - this is the graph info:
Colour 0 = reference


****************************************
  SUMMARY:
  Colour: MeanReadLen     TotalSeq
0       100     200
1       90      73500000000
****************************************
*/

  GraphInfo ginfo;
  graph_info_initialise(&ginfo);
  graph_info_set_seq(&ginfo, 0, 200);
  graph_info_set_seq(&ginfo, 1, 73500000000 );
  graph_info_set_mean_readlen(&ginfo, 0, 100);
  graph_info_set_mean_readlen(&ginfo, 1, 90);

  for (j=2; j<NUMBER_OF_COLOURS; j++)
    { 
      graph_info_set_seq(&ginfo, j, 0 );
      graph_info_set_mean_readlen(&ginfo, j, 0);
    }

  GraphAndModelInfo model_info;
  double mu=0.8;
  double seq_err_rate_per_base=0.01;
  int ref_colour=0;
  int num_chroms=2; 
  long long genome_len = 3000000000;
  initialise_model_info(&model_info, &ginfo, genome_len, mu, seq_err_rate_per_base, 
			ref_colour, num_chroms, EachColourADiploidSampleExceptTheRefColour);

  var.one_allele       = branch1;
  var.one_allele_or    =branch1_o;
  var.len_one_allele   = 56;
  var.other_allele     = branch2;
  var.other_allele_or    =branch2_o;
  var.len_other_allele = 56;
  var.which = first;
  int little_width =100;
  int little_height = 10;
  int little_retries=20;
  LittleHashTable* little_db_graph = little_hash_table_new(little_height, little_width, little_retries, db_graph->kmer_size);

  AnnotatedPutativeVariant annovar;
  GenotypingWorkingPackage* gwp = alloc_genotyping_work_package(60,1000,NUMBER_OF_COLOURS, NUMBER_OF_COLOURS+1);
  if (gwp==NULL)
    {
      printf("Unable to malloc genotyping work package. Out of memory. Abort\n");
      exit(1);
    }

  initialise_putative_variant(&annovar, &model_info, &var, SimplePathDivergenceCaller, 
			      55, AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice, gwp, db_graph, little_db_graph );




  CU_ASSERT(annovar.genotype[1]==hom_other);



  //******************
  /// second test

  fp = fopen("../data/test/pop_graph/variations/complex_genotyping/pd_example2_just_alleles.fasta", "r");
  if (fp==NULL)
    {
      printf("Unable to open test file ../data/test/pop_graph/variations/complex_genotyping/pd_example2_just_alleles.fasta\n");
      exit(1);
    }
  num_kmers_br1 = align_next_read_to_graph_and_return_node_array(fp, max_read_length, branch1, branch1_o, false, &file_reader_fasta, 
								      seq, kmer_window, db_graph, 0);

  num_kmers_br2 = align_next_read_to_graph_and_return_node_array(fp, max_read_length, branch2, branch2_o, false, &file_reader_fasta, 
								      seq, kmer_window, db_graph, 0);


  fclose(fp);
  //second meaning second test
  int second_br1_covg_0[]={1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  for (j=0; j<56; j++)
    {
      branch1[j]->coverage[0]=second_br1_covg_0[j];
    }

  int second_br1_covg_1[]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6};

  for (j=0; j<56; j++)
    {
      branch1[j]->coverage[1]=second_br1_covg_1[j];
    }



  int second_br2_covg_0[]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  for (j=0; j<56; j++)
    {
      branch2[j]->coverage[0]=second_br2_covg_0[j];
    }
  int second_br2_covg_1[]= {10,10,11,11,11,12,12,12,11,11,11,11,11,10,11,11,11,11,10,10,10,10,10,10,10,10,10,10,9,8,7,6,5,5,5,5,6,6,5,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,6,6};

  for (j=0; j<56; j++)
    {
      branch2[j]->coverage[1]=second_br2_covg_1[j];
    }
  


  AnnotatedPutativeVariant annovar2;
  wipe_little_graph(little_db_graph);
  initialise_putative_variant(&annovar2, &model_info, &var, SimplePathDivergenceCaller, 
			      55, AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice, gwp, db_graph, little_db_graph );


  CU_ASSERT(annovar2.genotype[1]==hom_other);//1==colour








  free_genotyping_work_package(gwp);
  little_hash_table_free(&little_db_graph);

  free(branch1);
  free(branch2);
  free(branch1_o);
  free(branch2_o);




}
