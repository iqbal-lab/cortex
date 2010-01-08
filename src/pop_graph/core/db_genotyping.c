#include <db_genotyping.h>


void set_variant_branches_and_flanks(VariantBranchesAndFlanks* var, 
				     dBNode** flank5p,    Orientation* flank5p_or,    int len_flank5p,
				     dBNode** ref_allele, Orientation* ref_allele_or, int len_ref_allele, 
				     dBNode** alt_allele, Orientation* alt_allele_or, int len_alt_allele, 
				     dBNode** flank3p,    Orientation* flank3p_or,    int len_flank3p)
{
  var->flank5p       = flank5p;
  var->flank5p_or    = flank5p_or;
  var->len_flank5p   = len_flank5p;
  var->ref_allele    = ref_allele;
  var->ref_allele_or = ref_allele_or;
  var->len_ref_allele= len_ref_allele;
  var->alt_allele    = alt_allele;
  var->alt_allele_or = alt_allele_or;
  var->len_alt_allele= len_alt_allele;
  var->flank3p       =flank3p;
  var->flank3p_or    = flank3p_or;
  var->len_flank3p   = len_flank3p;

}
				     


//Utility function - only exported so I can test it.
//Given two branches (arrays of dBNode*'s), return (in arguments 5,6,7,8) two arrays for each branch. One gives, for each node in the branch, the 
//number of times that node is seen in that branch. The second gives, for each node in a branch, the number of times that node is seen in
//the OTHER branch.
// The argument only_count_nodes_with_edge_in_specified_colour_func allows you to specify if you care whether a node has an edge in some colour or not.
// e.g. if we are only interested in the blue subgraph, and we do not want to count nodes that are not in the blue subgraph, we set this to true,
// and use the following two arguments also:
//The last two arguments allow you to do this to subgraph of the de Bruijn graph defined by these two functions - eg the union of all colours,
// or just one colour. These arguments are ignored if only_count_nodes_with_edge_in_specified_colour_func==false
void get_node_multiplicities(dBNode** branch1, int len_branch1, dBNode** branch2, int len_branch2, 
			     int** br1_multiplicity_in_br1, int** br2_multiplicity_in_br2,
			     int** br1_multiplicity_in_br2, int** br2_multiplicity_in_br1,
			     boolean only_count_nodes_with_edge_in_specified_colour_func,
			     Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*) )
{

  //we will reuse this;
  void get_mult(dBNode** br_src, int len_br_src, dBNode** br_target, int len_br_target, int** mult_array)
  {
    int i,j;
    int count_occurrences=0; //will be number of things we have put in this array

    for (i=0; i<len_br_src; i++)
      {
	*mult_array[i]=0;
      }

    for (i=0; i<len_br_src ; i++)
      {
	count_occurrences=0;
	
	for (j=0 ; j<len_br_target; j++)
	  {
	    //if i-th and j-th elements are the same, AND they exist in the colour (or function of colours) we are interested in
	    if (db_node_addr_cmp(&br_src[i], &br_target[j])==0 )
	      {
		if ( (only_count_nodes_with_edge_in_specified_colour_func==true) && 
		     (!db_node_is_this_node_in_subgraph_defined_by_func_of_colours(br_src[i], get_colour)) )
		  {
		    //does not count if node does not exist in the specified subgraph
		  }
		else
		  {
		    count_occurrences++;
		  }
	      }
	  }
	
	*mult_array[i]=count_occurrences;
      }
  }
  get_mult(branch1, len_branch1, branch1, len_branch1,  br1_multiplicity_in_br1);
  get_mult(branch2, len_branch2, branch2, len_branch2,  br2_multiplicity_in_br2);
  get_mult(branch1, len_branch1, branch2, len_branch2,  br1_multiplicity_in_br2);
  get_mult(branch2, len_branch2, branch1, len_branch1,  br2_multiplicity_in_br1);
  
}



void genotype_all_variants_in_fff_file(char* filename, Variant_File_Format output_format, char* out_filename, int max_branch_len,
				       void (*genotype_model) (VariantBranchesAndFlanks*, dBGraph*, int, Variant_File_Format, FILE*, 
							       int, int, int, double), 
				       dBGraph* db_graph, int colour1 , int colour_ref, int population_size_colour1, 
				       double avg_depth_of_covg_pop1, int avg_read_len)
{
  //call read_next_variant in a loop, and pass each one into the genotype_model function. The model function will print in the required format
  FILE* fptr = fopen(filename, "r");
  if (fptr==NULL)
    {
      printf("Unable to open %s\n", filename);
      exit(1);
    }
  FILE* output_fptr = fopen(out_filename, "w");
  if (fptr==NULL)
    {
      printf("Unable to open %s\n", out_filename);
      exit(1);
    }

  
  dBNode* flank5p[max_branch_len];//assume max_brnahc_len > max flank length
  dBNode* ref_allele[max_branch_len];
  dBNode* alt_allele[max_branch_len];
  dBNode* flank3p[max_branch_len];
  Orientation flank5p_or[max_branch_len];
  Orientation ref_allele_or[max_branch_len];
  Orientation alt_allele_or[max_branch_len];
  Orientation flank3p_or[max_branch_len];
  int len_flank5p;
  int len_ref_allele;
  int len_alt_allele;
  int len_flank3p;
  VariantBranchesAndFlanks var;

  while (read_next_variant_from_full_flank_file(fptr, max_branch_len, 
						flank5p, flank5p_or, &len_flank5p,
						ref_allele, ref_allele_or, &len_ref_allele,
						alt_allele, alt_allele_or, &len_alt_allele,
						flank3p, flank3p_or, &len_flank3p,
						db_graph, colour1))
    {

      set_variant_branches_and_flanks(&var, 
				      flank5p, flank5p_or, len_flank5p,
				      ref_allele, ref_allele_or, len_ref_allele,
				      alt_allele, alt_allele_or, len_alt_allele,
				      flank3p, flank3p_or, len_flank3p
				      );


      genotype_model(&var, db_graph, avg_read_len, glf, output_fptr, colour1, colour_ref, population_size_colour1, avg_depth_of_covg_pop1 );
    }
}





// average depth of covg is diploid coverage - we will sort out covg per allele internally to this function
void genotype_model_simple(VariantBranchesAndFlanks* var, dBGraph* db_graph, int avg_read_len,
			   Variant_File_Format output_format, FILE* output_fptr, 
			   int colour1 , int colour_ref, int population_size_colour1, double avg_depth_of_covg_pop1)
{

  int ref_multiplicity_in_ref[var->len_ref_allele];
  int ref_multiplicity_in_alt[var->len_ref_allele];
  int alt_multiplicity_in_alt[var->len_alt_allele];
  int alt_multiplicity_in_ref[var->len_alt_allele];
  int* ref_multiplicity_in_ref_ptr[var->len_ref_allele];
  int* ref_multiplicity_in_alt_ptr[var->len_ref_allele];
  int* alt_multiplicity_in_alt_ptr[var->len_alt_allele];
  int* alt_multiplicity_in_ref_ptr[var->len_alt_allele];
  int ref_normalised_covg[var->len_ref_allele];
  int alt_normalised_covg[var->len_alt_allele];
  int* ref_normalised_covg_ptr[var->len_ref_allele];
  int* alt_normalised_covg_ptr[var->len_alt_allele];
  
  int i;
  for (i=0; i<var->len_ref_allele; i++)
    {
      ref_multiplicity_in_ref[i]=0;
      ref_multiplicity_in_alt[i]=0;
      ref_multiplicity_in_ref_ptr[i]=&ref_multiplicity_in_ref[i];
      ref_multiplicity_in_alt_ptr[i]=&ref_multiplicity_in_alt[i];
      ref_normalised_covg[i]=0;
      ref_normalised_covg_ptr[i] = &ref_normalised_covg[i];
    }
  for (i=0; i<var->len_alt_allele; i++)
    {
      alt_multiplicity_in_ref[i]=0;
      alt_multiplicity_in_alt[i]=0;
      alt_multiplicity_in_ref_ptr[i]=&alt_multiplicity_in_ref[i];
      alt_multiplicity_in_alt_ptr[i]=&alt_multiplicity_in_alt[i];
      alt_normalised_covg[i]=0;
      alt_normalised_covg_ptr[i]=&alt_normalised_covg[i];
    }

  Edges get_colour1(const dBNode* node)
  {
    return  get_edge_copy(*node, individual_edge_array, colour1);
  }
  int get_covg1(const dBNode* node)
  {
    return node->coverage[colour1];
  }
  
  get_node_multiplicities(var->ref_allele, var->len_ref_allele, var->alt_allele, var->len_alt_allele, 
			  ref_multiplicity_in_ref_ptr,  alt_multiplicity_in_alt_ptr,
			  ref_multiplicity_in_alt_ptr,  alt_multiplicity_in_ref_ptr,
			  true,
			  &get_colour1, &get_covg1 );
  
  get_normalised_coverage(population_size_colour1, avg_depth_of_covg_pop1, avg_read_len, 
			  ref_normalised_covg, alt_normalised_covg,
			  ref_multiplicity_in_ref, alt_multiplicity_in_ref,
			  var, db_graph, colour1, colour_ref);
  
  double prob_data_under_alt_ref(int* mle_allele_freq) 
  {
    //will need to work through each branch, looking for maximal subcontigs that have multiplicity 1
    //(multiplicity refers to number of times seen on one branch or other)

    int i;
    double likelihoods[2*population_size_colour1];//likeihoods at different allele frequencies of ALT
    double max_likelihood=LARGE_NEGATIVE_NUMBER;
    int mle_af=-1;

    for (i=0; i<2*population_size_colour1; i++)
      {
	likelihoods[i]=calc_log_likelihood_of_data_seen_on_one_allele(ref_normalised_covg, ref_multiplicity_in_ref, ref_multiplicity_in_alt, ref_normalised_covg,
								      var->len_ref_allele, 2*population_size_colour1-i, avg_depth_of_covg_pop1/2, 
								      avg_read_len, db_graph->kmer_size)
	  + calc_log_likelihood_of_data_seen_on_one_allele(alt_normalised_covg, alt_multiplicity_in_alt, alt_multiplicity_in_ref, alt_normalised_covg,
									   var->len_alt_allele, i, avg_depth_of_covg_pop1/2, 
									   avg_read_len, db_graph->kmer_size);
	if (likelihoods[i]>mle_af)
	  {
	    max_likelihood = likelihoods[i];
	    mle_af = i;
	  }
      }
    *mle_allele_freq = mle_af;
    return max_likelihood;;

  }

  double prob_data_under_alt_alt(int* mle_allele_freq) 
  {
    //will need to work through each branch, looking for maximal subcontigs that have multiplicity 1
    //(multiplicity refers to number of times seen on one branch or other)

    int i;
    double likelihoods[2*population_size_colour1];//likeihoods at different allele frequencies of ALT
    double max_likelihood=LARGE_NEGATIVE_NUMBER;
    int mle_af=-1;

    for (i=0; i<2*population_size_colour1; i++)
      {
	likelihoods[i]= calc_log_likelihood_of_data_seen_on_one_allele(alt_normalised_covg, alt_multiplicity_in_alt, alt_multiplicity_in_ref, alt_normalised_covg,
								       var->len_alt_allele, i, avg_depth_of_covg_pop1, //note depth here is full diploid depth
								       avg_read_len, db_graph->kmer_size);
	if (likelihoods[i]>mle_af)
	  {
	    max_likelihood = likelihoods[i];
	    mle_af = i;
	  }
      }
    *mle_allele_freq = mle_af;
    return max_likelihood;;

  }

  int max_lik_af_AR;
  int max_lik_af_AA;

  fprintf(output_fptr, "%f\t%d\t%f\t%d\n", prob_data_under_alt_ref(&max_lik_af_AR), max_lik_af_AR, prob_data_under_alt_alt(&max_lik_af_AA), max_lik_af_AA);
}

