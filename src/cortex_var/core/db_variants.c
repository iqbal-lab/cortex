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

#include <math.h>
#include <db_variants.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <maths.h>
#include <stdlib.h>

void set_variant_branches_and_flanks(VariantBranchesAndFlanks* var, 
				     dBNode** flank5p,    Orientation* flank5p_or,    int len_flank5p,
				     dBNode** one_allele, Orientation* one_allele_or, int len_one_allele, 
				     dBNode** other_allele, Orientation* other_allele_or, int len_other_allele, 
				     dBNode** flank3p,    Orientation* flank3p_or,    int len_flank3p,
				     WhichAlleleIsRef which)
{
  var->flank5p       = flank5p;
  var->flank5p_or    = flank5p_or;
  var->len_flank5p   = len_flank5p;
  var->one_allele    = one_allele;
  var->one_allele_or = one_allele_or;
  var->len_one_allele= len_one_allele;
  var->other_allele    = other_allele;
  var->other_allele_or = other_allele_or;
  var->len_other_allele= len_other_allele;
  var->flank3p       =flank3p;
  var->flank3p_or    = flank3p_or;
  var->len_flank3p   = len_flank3p;
  var->which         = which;
}


/*
void print_both_alleles(VariantBranchesAndFlanks* var)
{
  
}
*/


void action_set_flanks_and_branches_to_be_ignored(VariantBranchesAndFlanks* var)
{
  int i;
  for (i=0; i<var->len_one_allele; i++)
    {
      db_node_action_set_status_ignore_this_node((var->one_allele)[i]);
    }
  for (i=0; i<var->len_other_allele; i++)
    {
      db_node_action_set_status_ignore_this_node((var->other_allele)[i]);
    }
  for (i=0; i<var->len_flank5p; i++)
    {
      db_node_action_set_status_ignore_this_node((var->flank5p)[i]);
    }
  for (i=0; i<var->len_flank3p; i++)
    {
      db_node_action_set_status_ignore_this_node((var->flank3p)[i]);
    }
}


void db_variant_action_do_nothing(VariantBranchesAndFlanks* var)
{
  return;
}


boolean  db_variant_precisely_one_allele_is_in_given_func_of_colours(VariantBranchesAndFlanks* var, Edges (*get_colour)(const dBNode*), dBGraph* db_graph, WhichAllele* which)

{
  
  if ( (does_this_path_exist_in_this_colour(var->one_allele, var->one_allele_or, var->len_one_allele, get_colour, db_graph)==true)
       &&
       (does_this_path_exist_in_this_colour(var->other_allele, var->other_allele_or, var->len_other_allele, get_colour,db_graph)==false)
       )
    {
      *which = allele_one;
      return true;
    }
  else if ( (does_this_path_exist_in_this_colour(var->one_allele, var->one_allele_or, var->len_one_allele, get_colour, db_graph)==false)
       &&
       (does_this_path_exist_in_this_colour(var->other_allele, var->other_allele_or, var->len_other_allele, get_colour,db_graph)==true)
       )
    {
      *which = allele_other;
      return true;
    } 
  else
    {
      return false;
    }
       
}

//very simplistic - returns hom if one allele is there in given colour, het if both are there, and absetnt if neither
zygosity db_variant_get_zygosity_in_given_func_of_colours(VariantBranchesAndFlanks* var, Edges (*get_colour)(const dBNode*), dBGraph* db_graph)
{
  if (does_this_path_exist_in_this_colour(var->one_allele, var->one_allele_or, var->len_one_allele, get_colour, db_graph)==true)       
    {
      if (does_this_path_exist_in_this_colour(var->other_allele, var->other_allele_or, var->len_other_allele, get_colour,db_graph)==true)
	{
	  return het;
	}
      else
	{
	  return hom_one;
	}
    }
  else if (does_this_path_exist_in_this_colour(var->other_allele, var->other_allele_or, var->len_other_allele, get_colour,db_graph)==true)
    {
      return hom_other;
    }
  else
    {
      return absent;
    }
  
}


double get_log_bayes_factor_comparing_genotypes_at_bubble_call(zygosity genotype1, zygosity genotype2, VariantBranchesAndFlanks* var, 
							      double seq_error_rate_per_base, double sequencing_depth_of_coverage, int read_length, int colour)
{
  boolean too_short = false;
  int initial_covg_plus_upward_jumps_branch1 = count_reads_on_allele_in_specific_colour(var->one_allele, var->len_one_allele, colour, &too_short);
  int initial_covg_plus_upward_jumps_branch2 = count_reads_on_allele_in_specific_colour(var->other_allele, var->len_other_allele, colour, &too_short);
  if (too_short==true)
    {
      return 0;
    }

  double theta_one = ((double)(sequencing_depth_of_coverage * var->len_one_allele))/( (double) read_length );
  double theta_other = ((double)(sequencing_depth_of_coverage * var->len_other_allele))/( (double) read_length );

  double log_prob_data_given_genotype1 = get_log_likelihood_of_genotype_under_poisson_model_for_covg_on_variant_called_by_bubblecaller(genotype1, seq_error_rate_per_base, 
																       initial_covg_plus_upward_jumps_branch1, 
																       initial_covg_plus_upward_jumps_branch2, theta_one, theta_other, var);
  double log_prob_data_given_genotype2 = get_log_likelihood_of_genotype_under_poisson_model_for_covg_on_variant_called_by_bubblecaller(genotype2, seq_error_rate_per_base, 
																       initial_covg_plus_upward_jumps_branch1, 
																       initial_covg_plus_upward_jumps_branch2, theta_one, theta_other, var);

  return log_prob_data_given_genotype1 - log_prob_data_given_genotype2;
}


//assuming a pair of branches really do make up a variant, calculate the log likelihood of a genotype
//under the model described in our paper (eg used for HLA)
//theta here is as used in the paper: (D/R) * length of branch/allele. NOT the same theta as seen in model_selection.c
//assumes called by BubbleCaller, so no overlaps between alleles.
double get_log_likelihood_of_genotype_under_poisson_model_for_covg_on_variant_called_by_bubblecaller(zygosity genotype, double error_rate_per_base, int covg_branch_1, int covg_branch_2, 
												     double theta_one, double theta_other, VariantBranchesAndFlanks* v)
{

  if (genotype==hom_one)
    {
      return 2* (covg_branch_1*log(theta_one/2) - theta_one/2 - log_factorial(covg_branch_1) + covg_branch_2 *log(error_rate_per_base) );
    }
  else if (genotype==hom_other)
    {
      return 2* (covg_branch_2*log(theta_other/2) - theta_other/2 - log_factorial(covg_branch_2) + covg_branch_1 *log(error_rate_per_base)  );
    }
  else if (genotype==het)
    {
      return (covg_branch_1*log(theta_one/2) - theta_one/2 - log_factorial(covg_branch_1) )
        + (covg_branch_2*log(theta_other/2) - theta_other/2 - log_factorial(covg_branch_2) );
    }
  else
    {
      printf("Programming error. called get_log_likelihood_of_genotype_under_poisson_model_for_covg_on_variant_called_by_bubblecaller with bad genotype");
      exit(1);
    }
}










//returns true if was able to initialise
// if you have no GraphInfo, enter NULL
//if you do not know what sequencing_error_per_base is, enter -1, will use default of 0.01
//if you do not knwo what genome length is, enter -1, will use default of 3 billion (human)
boolean initialise_putative_variant(AnnotatedPutativeVariant* annovar, VariantBranchesAndFlanks* var, DiscoveryMethod caller, 
				    GraphInfo* ginfo, double seq_error_rate_per_base, long long genome_length)
{
  annovar->caller = caller;
  annovar->var=var;

  boolean flag1=false;
  boolean flag2=false;

  flag1=get_num_effective_reads_on_branch(annovar->br1_covg, var->one_allele, var->len_one_allele);
  flag2=get_num_effective_reads_on_branch(annovar->br2_covg, var->other_allele, var->len_other_allele);
  
  if (var->len_one_allele < var->len_other_allele)
    {
      annovar->len_start=var->len_one_allele;
    }
  else
    {
      annovar->len_start=var->len_other_allele;
    }
  if (annovar->len_start<2)
    {
      annovar->too_short = true;
    }
  else
    {
      annovar->too_short = false;
    }

  if (    (( (flag1==false) || (flag2==false) ) && (annovar->too_short==true))
	  ||
	  (( (flag1==true) || (flag2==true) ) && (annovar->too_short==false))
	  )
    {
      printf("Unexpected inconsistency in initialise_putative_variant - coding error\n");
      exit(1);
    }

 get_num_effective_reads_on_branch(annovar->theta1, var->one_allele, annovar->len_start);
 get_num_effective_reads_on_branch(annovar->theta2, var->other_allele, annovar->len_start);
  
  int i;
  annovar->BigTheta = 0;
  annovar->BigThetaStart = 0;
  for (i=0; i<var->len_one_allele; i++)
    {
      annovar->BigTheta += annovar->br1_covg[i];
    }
  for (i=0; i<var->len_other_allele; i++)
    {
      annovar->BigTheta += annovar->br2_covg[i];
    }

  for (i=0; i<annovar->len_start; i++)
    {
      annovar->BigThetaStart += annovar->theta1[i] + annovar->theta2[i];
    }

  //determine genotype (under assumption/model that this is a variant) for each colour
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      double sequencing_depth_of_coverage=0;
      if (seq_error_rate_per_base==-1)
	{
	  seq_error_rate_per_base=0.01;//default
	}
      if (genome_length==-1)
	{
	  genome_length=3000000000;//default
	}
      int mean_read_len=100;
      if (ginfo !=NULL)
	{
	  sequencing_depth_of_coverage=(double) ginfo->total_sequence[i]/genome_length;
	  mean_read_len = ginfo->mean_read_length[i];
	}
      if (sequencing_depth_of_coverage==0)
	{
	  annovar->genotype[i]=absent;
	}
      else
	{
	  double log_bf_hom1_over_het= get_log_bayes_factor_comparing_genotypes_at_bubble_call(hom_one, het,  var, seq_error_rate_per_base, sequencing_depth_of_coverage,
											     mean_read_len,i);
	  double log_bf_hom1_over_hom2= get_log_bayes_factor_comparing_genotypes_at_bubble_call(hom_one, hom_other,  var, seq_error_rate_per_base, sequencing_depth_of_coverage,
											      mean_read_len,i);
	  double log_bf_het_over_hom2= get_log_bayes_factor_comparing_genotypes_at_bubble_call(het, hom_other,  var, seq_error_rate_per_base, sequencing_depth_of_coverage,
											     mean_read_len,i);
	  
	  if ( ( log_bf_hom1_over_het >0 ) && (log_bf_hom1_over_hom2>0) )
	    {
	      annovar->genotype[i]=hom_one;
	    }
	  else if ( ( log_bf_hom1_over_het <0 ) && (log_bf_het_over_hom2>0) )
	    {
	      annovar->genotype[i]=het;
	    }
	  else if ( ( log_bf_hom1_over_het <0 ) && (log_bf_het_over_hom2<0) )
	    {
	      annovar->genotype[i]=hom_other;
	    }
	  else
	    {
	      printf("Coding error - get incompatible combination of bates factors: %f, %f, %f\n", log_bf_hom1_over_het, log_bf_hom1_over_hom2, log_bf_hom1_over_het);
	      exit(1);
	    }
	}
    }

  return true;
}

//first argument is an array of length NUMBER_OF_COLOURS, into which results go.
//If you want the number of reads on the entire branch, enter the length of that branch in arg3 (eg var->len_one-allele)
//Sometimes we want to take just the start of the branch (if one branch is longer than the other, we may just take the length of the shorter one)
//and so you enter that in arg3 in that case
//note these are effective reads, as counting covg in the de Bruijn graph
//returns false if branch is too short (1 or 2 nodes) to do this
boolean get_num_effective_reads_on_branch(int* array, dBNode** allele, int how_many_nodes)
{
  int i;
  boolean too_short=false;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      array[i] = count_reads_on_allele_in_specific_colour(allele, how_many_nodes, i, &too_short);
      if (too_short==true)
	{
	  return too_short;
	}
    }
  return too_short;
}




//does not count covg on first or last nodes, as they are bifurcation nodes
//if length==0 or 1  returns -1.
int count_reads_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short)
{

  if ( (len==0) || (len==1) )
    {
      *too_short=true;
      return -1;
    }

  //note start at node 1, avoid first node
  int total= db_node_get_coverage(allele[1], individual_edge_array, colour);

  int i;

  //note we do not go as far as the final node, which is where the two branches rejoin
  for (i=2; i<len-1; i++)
    {
      if (db_node_get_coverage(allele[i], individual_edge_array, i) 
	  > db_node_get_coverage(allele[i-1], individual_edge_array, i-1) )
        {
          total++;
        }
    }

  return total;
}

