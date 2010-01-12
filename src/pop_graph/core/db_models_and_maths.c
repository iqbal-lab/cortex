#include <db_models_and_maths.h>

const double LARGE_NEGATIVE_NUMBER=-1000000.0;


int factorial( int n )
{
  if ( n <= 1 )
    return 1;
  else
    return  n * factorial( n-1 );
}

int max(int a, int b, int c, int d)
{
  if ( (a>=b) && (a>=c) && (a>=d) )
    {
      return a;
    }
  else if ( (b>=a) && (b>=c) && (b>=d) )
    {
      return b;
    }
  else if ( (c>=a) && (c>=b) && (c>=d) )
    {
      return c;
    }
  else
    {
      return d;
    }
}



// Given a sequence of coverage values for a set of contiguous nodes,
// where we have factored out effects due to rest of genome or repeat/multiplicity
//calculate probability, under a Poisson approx to a binomial model, 
// that we see those coverages
double get_probability_of_observed_covgs(int* coverages, int array_len, double lambda)
{
  //Get positions where coverage jumps up.
  int jumps[array_len];
  int i;
  int num_jumps=0;
  jumps[num_jumps]=0;

  for (i=1; i<array_len; i++)
    {
      if (coverages[i]>coverages[i-1])
	{
	  num_jumps++;
	  jumps[num_jumps]=i;//note I am putting the k-th jump at index k
	}
    }

  //Get probability that a Poisson(lambda) variable takes value coverages[0]
  double prob_start = exp(-lambda)*pow(lambda, (double)coverages[0] )/factorial(coverages[0]);

  //Get probability that a series of Exponential(lambda) variables could give the sequence of jumps
  double prob_jumps = 1;
  for (i=1; i<=num_jumps; i++)
    {
      prob_jumps=prob_jumps*( exp(jumps[i]-jumps[i-1]-1) - exp(jumps[i]-jumps[i-1]) );
    }

  return prob_start*prob_jumps;
}

// Given a sequence of coverage values for a set of contiguous nodes,
// where we have factored out effects due to rest of genome or repeat/multiplicity
//calculate log probability, under a Poisson approx to a binomial model, 
// that we see those coverages
double get_log_probability_of_observed_covgs(int* coverages, int array_len, double lambda)
{
  //Get positions where coverage jumps up.
  int jumps[array_len];
  int i;
  int num_jumps=0;
  jumps[num_jumps]=0;

  for (i=1; i<array_len; i++)
    {
      if (coverages[i]>coverages[i-1])
	{
	  num_jumps++;
	  jumps[num_jumps]=i;//note I am putting the k-th jump at index k
	}
    }

  //Get log probability that a Poisson(lambda) variable takes value coverages[0]
  double log_prob_start = coverages[0]*log(lambda) -lambda -log((double)factorial(coverages[0])); 

  //Get log probability that a series of Exponential(lambda) variables could give the sequence of jumps
  double log_prob_jumps = 1;
  for (i=1; i<=num_jumps; i++)
    {
      log_prob_jumps=log_prob_jumps+ log( exp( (double) (jumps[i]-jumps[i-1]-1) ) - exp( (double) (jumps[i]-jumps[i-1]))  );
    }

  return log_prob_start+log_prob_jumps;
}


//normalise away the coverage due to the rest of the genome. This requires you to know which allele is the ref allele
//If you do not know that, use the get_rough_normalised_coverage which does a rough approximation
void get_normalised_coverage(int population_size, int avg_depth_of_covg, int avg_read_len, 
			    int* ref_normalised_covg, int*alt_normalised_covg,
			    int* ref_multiplicity_in_ref, int* alt_multiplicity_in_ref,
			    VariantBranchesAndFlanks* var, dBGraph* db_graph, int colour, int colour_ref)
{
  
  double diploid_lambda  = population_size*avg_depth_of_covg*(1-(db_graph->kmer_size)/avg_read_len);
  
  if (var->which != first)
    {
      printf("Cannot call get_normalised_coverage unless you can guarantee the ref allele is the first in var, and this is reflected in all your multiplicity arrays etc");
      exit(1);
    }

  //get normalised coverages
  int i;
  for (i=0; i<var->len_one_allele; i++)
    {
      //normalised covg = covg of data in colour1 - diploid_lambda(number of times seem in ref genome - num times seen in ref allele)
      ref_normalised_covg[i] = db_node_get_coverage((var->one_allele)[i], individual_edge_array, colour) - 
	diploid_lambda*( db_node_get_coverage((var->one_allele)[i], individual_edge_array, colour_ref) - ref_multiplicity_in_ref[i]);
    }
  
  for (i=0; i<var->len_other_allele; i++)
    {
      //normalised covg = covg of data in colour1 - diploid_lambda(number of times seem in ref genome - num times seen in ref allele)
      alt_normalised_covg[i] = db_node_get_coverage((var->other_allele)[i], individual_edge_array, colour) - 
	diploid_lambda*( db_node_get_coverage((var->one_allele)[i], individual_edge_array, colour_ref) - alt_multiplicity_in_ref[i]);
    }
}




void get_rough_normalised_coverage(int population_size, int avg_depth_of_covg, int avg_read_len, 
				   int* ref_normalised_covg, int*alt_normalised_covg,
				   int* ref_multiplicity_in_ref, int* alt_multiplicity_in_alt,
				   VariantBranchesAndFlanks* var, dBGraph* db_graph, int colour, int colour_ref)
{

  if (var->which != unknown)
    {
      printf("Why are you calling get_rough_normalised_coverage when you KNOW which allele is the ref?\n");
      exit(1);
    }

  int i;
  for (i=0; i<var->len_one_allele; i++)
    {
      int multiplicity_in_ref_genome = db_node_get_coverage((var->one_allele)[i], individual_edge_array, colour_ref);
      
      if (multiplicity_in_ref_genome>1)
	{
	  //normalised covg = covg of data in colour1 - pop_size*(multiplicity in ref - multiplicity in this allele)*avg_covg_per_person*scaling factor because of kmer_length&read length
	  ref_normalised_covg[i] = db_node_get_coverage((var->one_allele)[i], individual_edge_array, colour) - 
	    (multiplicity_in_ref_genome - ref_multiplicity_in_ref[i])*population_size*avg_depth_of_covg*(1-(db_graph->kmer_size)/avg_read_len);
	}
      else
	{
	  ref_normalised_covg[i]=db_node_get_coverage((var->one_allele)[i], individual_edge_array, colour);
	}
    }
  
  for (i=0; i<var->len_other_allele; i++)
    {
      int multiplicity_in_ref_genome = db_node_get_coverage((var->other_allele)[i], individual_edge_array, colour_ref);

      if (multiplicity_in_ref_genome>1)
	{
	  alt_normalised_covg[i] = db_node_get_coverage((var->other_allele)[i], individual_edge_array, colour) - 
	    (multiplicity_in_ref_genome - alt_multiplicity_in_alt[i])*population_size*avg_depth_of_covg*(1-(db_graph->kmer_size)/avg_read_len);
	}
      else
	{
          alt_normalised_covg[i]=db_node_get_coverage((var->other_allele)[i], individual_edge_array, colour);
        }

    }

}


//allele freq is int between 0 and 2*population_size
//This calculates, conditional on the allele frequency, the log likelihood of seeing the data we see on that allele
//Uses a Poisson model, but also breaks the allele into intervals consisting of nodes that are unique on that allele
//Then calculates the likelihood for the union of those, and multiplies by the proportion of the allele that lies in these intervals
// i.e. if there are 5 nodes that are repeated twice within the allele, and 200 that are unique, then it will do the calculation for the 200,
// and multiply by 200/205. (and then return the log of this)
// Note, if you have a bunch of people sequenced to 5x depth each, the avg_depth_of_covg_per_haploid is 2.5. 
double calc_log_likelihood_of_data_seen_on_one_allele_excluding_nodes_on_both_alleles(int* array_of_covgs, int* array_of_mult_wrt_self, int* array_of_mult_wrt_other, int* normalised_covg, 
										      int len_arrays, int allele_freq, double avg_depth_of_covg_per_haploid, int avg_read_len, short kmer_size)
{
    double haploid_lambda  = allele_freq*avg_depth_of_covg_per_haploid*(1-kmer_size/avg_read_len);
    
    int array[len_arrays];
    int array_count=0;
    int k=0;
    int num_nodes_with_mult_above1=0;
    double log_prob_data=0;

    while (k<len_arrays)
      {
	if ((array_of_mult_wrt_self[k]>1) || (array_of_mult_wrt_other[k]>0) )//this node occurs >1 time in this allele, or >0 times in the other allele
	  {
	    if (array_count>0)
	      {
		log_prob_data += get_log_probability_of_observed_covgs(array, array_count, haploid_lambda); 
	      }
	    k++;
	    array_count=0;
	    num_nodes_with_mult_above1++;
	  }
	else
	  {
	    //still in the same contiguous chunk of nodes with multiplicity 1 in the ref allele
	    array[array_count]=normalised_covg[k];
	    k++;
	    array_count++;
	  }
      }
    double proportion_of_ref_allele_with_multiplicity1 = (len_arrays - num_nodes_with_mult_above1)/len_arrays;
    return log_prob_data+ log(proportion_of_ref_allele_with_multiplicity1);
    
    
}



//allele freq is int between 0 and 2*population_size
//This calculates, conditional on the allele frequency, the log likelihood of seeing the data we see on that allele
//Uses a Poisson model, but also breaks the allele into intervals consisting of nodes that are unique on that allele
//Then calculates the likelihood for the union of those, and multiplies by the proportion of the allele that lies in these intervals
// i.e. if there are 5 nodes that are repeated twice within the allele, and 200 that are unique, then it will do the calculation for the 200,
// and multiply by 200/205. (and then return the log of this)
// Note, if you have a bunch of people sequenced to 5x depth each, the avg_depth_of_covg_per_haploid is 2.5. 
double calc_log_likelihood_of_data_seen_on_one_allele(int* array_of_covgs, int* array_of_mult_wrt_self,  int* normalised_covg, 
						      int len_arrays, int allele_freq, double avg_depth_of_covg_per_haploid, int avg_read_len, short kmer_size)
{
    double haploid_lambda  = allele_freq*avg_depth_of_covg_per_haploid*(1-kmer_size/avg_read_len);
    
    int array[len_arrays];
    int array_count=0;
    int k=0;
    int num_nodes_with_mult_above1=0;
    double log_prob_data=0;

    while (k<len_arrays)
      {
	if (array_of_mult_wrt_self[k]>1)//this node occurs >1 time in this allele
	  {
	    if (array_count>0)
	      {
		log_prob_data += get_log_probability_of_observed_covgs(array, array_count, haploid_lambda); 
	      }
	    k++;
	    array_count=0;
	    num_nodes_with_mult_above1++;
	  }
	else
	  {
	    //still in the same contiguous chunk of nodes with multiplicity 1 in the ref allele
	    array[array_count]=normalised_covg[k];
	    k++;
	    array_count++;
	  }
      }
    double proportion_of_ref_allele_with_multiplicity1 = (len_arrays - num_nodes_with_mult_above1)/len_arrays;
    return log_prob_data+ log(proportion_of_ref_allele_with_multiplicity1);
    
    
}


