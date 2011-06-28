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
#include <gsl_sf_gamma.h>


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


void get_all_genotype_log_likelihoods_at_bubble_call_for_one_colour(AnnotatedPutativeVariant* annovar, double seq_error_rate_per_base, double sequencing_depth_of_coverage, int read_length, int colour)
{
  boolean too_short = false;
  int initial_covg_plus_upward_jumps_branch1 = count_reads_on_allele_in_specific_colour(annovar->var->one_allele, annovar->var->len_one_allele, colour, &too_short);
  int initial_covg_plus_upward_jumps_branch2 = count_reads_on_allele_in_specific_colour(annovar->var->other_allele, annovar->var->len_other_allele, colour, &too_short);

  if (too_short==true)
    {
      annovar->too_short=true;
      int j;
      //set all the log likelihoods to zero.
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  initialise_genotype_log_likelihoods(&(annovar->gen_log_lh[j]));
	}
      return ;
    }

  double theta_one = ((double)(sequencing_depth_of_coverage * annovar->var->len_one_allele))/( (double) read_length );
  double theta_other = ((double)(sequencing_depth_of_coverage * annovar->var->len_other_allele))/( (double) read_length );


  annovar->gen_log_lh[colour].log_lh[hom_one]   = 
    get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, seq_error_rate_per_base, 
								     initial_covg_plus_upward_jumps_branch1, 
								     initial_covg_plus_upward_jumps_branch2, 
								     theta_one, theta_other, annovar->kmer);
  annovar->gen_log_lh[colour].log_lh[het]       = 
    get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(het, seq_error_rate_per_base, 
								     initial_covg_plus_upward_jumps_branch1, 
								     initial_covg_plus_upward_jumps_branch2, 
								     theta_one, theta_other, annovar->kmer);

  annovar->gen_log_lh[colour].log_lh[hom_other] = 
    get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_other, seq_error_rate_per_base, 
								     initial_covg_plus_upward_jumps_branch1, 
								     initial_covg_plus_upward_jumps_branch2, 
								     theta_one, theta_other, annovar->kmer);
 
//  printf("Log likelihood of data in colour %d under hom_one is %f\n", colour, annovar->gen_log_lh[colour].log_lh[hom_one]);
//  printf("Log likelihood of data in colour %d under het is %f\n", colour, annovar->gen_log_lh[colour].log_lh[het] );
//  printf("Log likelihood of dat ain colour %d under hom_other is %f\n", colour, annovar->gen_log_lh[colour].log_lh[hom_other] );


}





void get_all_haploid_genotype_log_likelihoods_at_bubble_call_for_one_colour(AnnotatedPutativeVariant* annovar, double seq_error_rate_per_base, double sequencing_depth_of_coverage, int read_length, int colour)
{
  boolean too_short = false;
  int initial_covg_plus_upward_jumps_branch1 = count_reads_on_allele_in_specific_colour(annovar->var->one_allele, annovar->var->len_one_allele, colour, &too_short);
  int initial_covg_plus_upward_jumps_branch2 = count_reads_on_allele_in_specific_colour(annovar->var->other_allele, annovar->var->len_other_allele, colour, &too_short);

  if (too_short==true)
    {
      annovar->too_short=true;
      int j;
      //set all the log likelihoods to zero.
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  initialise_genotype_log_likelihoods(&(annovar->gen_log_lh[j]));
	}
      return ;
    }

  double theta_one = ((double)(sequencing_depth_of_coverage * annovar->var->len_one_allele))/( (double) read_length );
  double theta_other = ((double)(sequencing_depth_of_coverage * annovar->var->len_other_allele))/( (double) read_length );


  annovar->gen_log_lh[colour].log_lh[hom_one]   = 
    get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_one, seq_error_rate_per_base, 
								     initial_covg_plus_upward_jumps_branch1, 
								     initial_covg_plus_upward_jumps_branch2, 
								     theta_one, theta_other, annovar->kmer);

  annovar->gen_log_lh[colour].log_lh[hom_other] = 
    get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(hom_other, seq_error_rate_per_base, 
								     initial_covg_plus_upward_jumps_branch1, 
								     initial_covg_plus_upward_jumps_branch2, 
								     theta_one, theta_other, annovar->kmer);

  //this is not allowed by the model
  annovar->gen_log_lh[colour].log_lh[het]       = annovar->gen_log_lh[colour].log_lh[hom_one] + annovar->gen_log_lh[colour].log_lh[hom_other] - 9999999999;



}


//assuming a pair of branches really do make up a variant, calculate the log likelihood of a genotype
//under the model described in our paper (eg used for HLA)
//theta here is as used in the paper: (D/R) * length of branch/allele. NOT the same theta as seen in model_selection.c
//assumes called by BubbleCaller, so no overlaps between alleles.
double get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(zygosity genotype, double error_rate_per_base, int covg_branch_1, int covg_branch_2, 
									double theta_one, double theta_other, int kmer)
{

  if (genotype==hom_one)
    {
      //Apply formula for likelihood in section 9.0 of Supp. Methods of paper; no unique segment, one shared segment
      return covg_branch_1*log(theta_one) - theta_one - gsl_sf_lnfact(covg_branch_1)  + covg_branch_2 *log(error_rate_per_base) ;
    }
  else if (genotype==hom_other)
    {
      //Apply formula for likelihood in section 9.0 of Supp. Methods of paper; no unique segment, one shared segment
      return covg_branch_2*log(theta_other) - theta_other - gsl_sf_lnfact(covg_branch_2) + covg_branch_1 *log(error_rate_per_base)  ;
    }
  else if (genotype==het)
    {
      //Apply formula for likelihood in section 9.0 of Supp. Methods of paper; no shared segment, TWO unique segments
      return (covg_branch_1*log(theta_one/2) - theta_one/2 - gsl_sf_lnfact(covg_branch_1) )
        + (covg_branch_2*log(theta_other/2) - theta_other/2 - gsl_sf_lnfact(covg_branch_2) );
    }
  else
    {
      printf("Programming error. called get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller with bad genotype");
      exit(1);
    }
}







void initialise_genotype_log_likelihoods(GenotypeLogLikelihoods* gl)
{
  int i;
  gl->log_lh[hom_one]=0;
  gl->log_lh[het]=0;
  gl->log_lh[hom_other]=0;
}


//returns true if was able to initialise
// if you have no GraphInfo, enter NULL
//if you do not know what sequencing_error_per_base is, enter -1, will use default of 0.01
//if you do not knwo what genome length is, enter -1, will use default of 3 billion (human)
// NOTE THIS GENOTYPES THE SITE
// if you enter ref_colour=-1, there is no ref. Otherwise, the ref colour will be ignored for calculations like theta, where
// we are aggregating covg.
boolean initialise_putative_variant(AnnotatedPutativeVariant* annovar, VariantBranchesAndFlanks* var, DiscoveryMethod caller, 
				    GraphInfo* ginfo, double seq_error_rate_per_base, long long genome_length, int kmer,
				    int ref_colour, ExperimentType expt)
{
  int number_individals=NUMBER_OF_COLOURS;
  if ( (ref_colour !=-1) && ( (ref_colour<0) || (ref_colour>=NUMBER_OF_COLOURS) ) )
    {
      printf("Called initialise_putative_variant with a ref_colour which is not -1, nor one of the colours this executable is compiled for. Coding error - call Zam\n");
      exit(1);
    }
  if (ref_colour != -1)
    {
      number_individals--;
    }

  annovar->kmer = kmer;
  annovar->caller = caller;
  annovar->var=var;

  boolean flag1=false;
  boolean flag2=false;

  flag1=get_num_effective_reads_on_branch(annovar->br1_covg, var->one_allele, var->len_one_allele);
  flag2=get_num_effective_reads_on_branch(annovar->br2_covg, var->other_allele, var->len_other_allele);

  if ( (flag1==true)||(flag2==true) )
    {
      annovar->too_short = true;
    }
  else
    {
      annovar->too_short = false;
    }
  
  if (var->len_one_allele < var->len_other_allele)
    {
      annovar->len_start=var->len_one_allele ;
    }
  else
    {
      annovar->len_start=var->len_other_allele;
    }

  if (annovar->too_short==true)
    {
      //deliberately not setting values to zero, for performance reasons.
      //this annovar should be discarded as soon as this function returns.
    }
  else//is a valid putative site
    {
      get_num_effective_reads_on_branch(annovar->theta1, var->one_allele, annovar->len_start);
      get_num_effective_reads_on_branch(annovar->theta2, var->other_allele, annovar->len_start);
      
      int i;
      annovar->BigTheta = 0;
      annovar->BigThetaStart = 0;
      
      for (i=0; i<var->len_one_allele; i++)
	{
	  if (i==ref_colour)
	    {
	      continue;
	    }
	  annovar->BigTheta += annovar->br1_covg[i];
	}
      
      for (i=0; i<var->len_other_allele; i++)
	{
	  annovar->BigTheta     += annovar->br2_covg[i];
	}
      
      
      for (i=0; i<annovar->len_start; i++)
	{
	  annovar->BigThetaStart += annovar->theta1[i] + annovar->theta2[i];
	}
      
      //determine genotype (under assumption/model that this is a variant) for each colour
      
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  initialise_genotype_log_likelihoods(&(annovar->gen_log_lh[i]));
	}
      
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
	      printf("ZAM test - I think this check can be removed - this code should never be hit\n");
	    }
	  int mean_read_len=100;
	  if (ginfo !=NULL)
	    {
	      sequencing_depth_of_coverage=(double) ginfo->total_sequence[i]/genome_length;
	      mean_read_len = ginfo->mean_read_length[i];
	    }
	  if ( (sequencing_depth_of_coverage==0)|| (annovar->too_short==true) )
	    {
	      annovar->genotype[i]=absent;
	    }
	  else if ( (expt==EachColourADiploidSample) || (expt==EachColourADiploidSampleExceptTheRefColour) )
	    {
	      
	      get_all_genotype_log_likelihoods_at_bubble_call_for_one_colour(annovar,  seq_error_rate_per_base, sequencing_depth_of_coverage,mean_read_len,i);
	      
	      if (annovar->gen_log_lh[i].log_lh[hom_one]>= annovar->gen_log_lh[i].log_lh[het])
		{
		  if (annovar->gen_log_lh[i].log_lh[hom_one]>=annovar->gen_log_lh[i].log_lh[hom_other])
		    {
		      annovar->genotype[i]=hom_one;
		    }
		  else
		    {
		      annovar->genotype[i]=hom_other;
		    }
		}
	      else if (annovar->gen_log_lh[i].log_lh[het]>=annovar->gen_log_lh[i].log_lh[hom_other])
		{
		  annovar->genotype[i]=het;
		}
	      else
		{
		  annovar->genotype[i]=hom_other;
		}
	      
	    }
	  else if ( (expt==EachColourAHaploidSample) || (expt==EachColourAHaploidSampleExceptTheRefColour) )
	    {
	      get_all_haploid_genotype_log_likelihoods_at_bubble_call_for_one_colour(annovar,  seq_error_rate_per_base, sequencing_depth_of_coverage,mean_read_len,i);
	      
	      
	      if (annovar->gen_log_lh[i].log_lh[hom_one]> annovar->gen_log_lh[i].log_lh[hom_other])
		{
		  annovar->genotype[i]=hom_one;
		}
	      else
		{
		  annovar->genotype[i]=hom_other;
		}

	    }
	  else
	    {
	      annovar->genotype[i]=hom_one;
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
//returns TRUE if branch is too short (1 or 2 nodes) to do this
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
      int jump = db_node_get_coverage(allele[i], individual_edge_array, colour) - db_node_get_coverage(allele[i-1], individual_edge_array, colour);

      //we add a little check to ensure that we ignore isolated nodes with higher covg - we count jumps only if they are signifiers of a new read arriving
      // and one node does not a read make
      int diff_between_next_and_prev=-1;
      if (i<len-2)
	{
	  diff_between_next_and_prev=    db_node_get_coverage(allele[i+1], individual_edge_array, colour)
	                               - db_node_get_coverage(allele[i-1], individual_edge_array, colour);
	}
      
      if ( (jump>0) && (diff_between_next_and_prev!=0) )
        {
	      total=total+jump;
        }
    }
  return total;
}


//does not count covg on first or last nodes, as they are bifurcation nodes
//if length==0 or 1  returns -1.
int count_reads_on_allele_in_specific_func_of_colours(dBNode** allele, int len, int (*sum_of_covgs_in_desired_colours)(const Element *), boolean* too_short)
{

  if ( (len==0) || (len==1) )
    {
      *too_short=true;
      return -1;
    }

  //note start at node 1, avoid first node
  int total= sum_of_covgs_in_desired_colours(allele[1]);

  int i;

  //note we do not go as far as the final node, which is where the two branches rejoin
  for (i=2; i<len-1; i++)
    {
      int jump = sum_of_covgs_in_desired_colours(allele[i]) - sum_of_covgs_in_desired_colours(allele[i-1]);

      //we add a little check to ensure that we ignore isolated nodes with higher covg - we count jumps only if they are signifiers of a new read arriving
      // and one node does not a read make
      int diff_between_next_and_prev=-1;
      if (i<len-2)
	{
	  diff_between_next_and_prev=    sum_of_covgs_in_desired_colours(allele[i+1])- sum_of_covgs_in_desired_colours(allele[i-1]);
	}
      
      if ( (jump>0) && (diff_between_next_and_prev!=0) )
        {
	      total=total+jump;
        }
    }
  return total;
}

