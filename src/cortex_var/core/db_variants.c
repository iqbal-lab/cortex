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

VariantBranchesAndFlanks* alloc_VariantBranchesAndFlanks_object(int len_5p, int len_br1, int len_br2, int len_3p, int kmer_size)
{
  VariantBranchesAndFlanks* var = (VariantBranchesAndFlanks*) malloc(sizeof(VariantBranchesAndFlanks));
  if (var==NULL)
    {
      return var;
    }
  else
    {      
      var->flank5p      = (dBNode**) malloc(sizeof(dBNode*)*(len_5p));
      if (var->flank5p==NULL) 
	{
	  return NULL;
	}

      var->one_allele   = (dBNode**) malloc(sizeof(dBNode*)*(len_br1));
      if (var->one_allele==NULL) 
	{
	  free(var->flank5p);
	  free(var);
	  return NULL;
	}

      var->other_allele = (dBNode**) malloc(sizeof(dBNode*)*(len_br2));
      if (var->other_allele==NULL) 
	{
	  free(var->flank5p);
	  free(var->one_allele);
	  free(var);
	  return NULL;
	}

      var->flank3p      = (dBNode**) malloc(sizeof(dBNode*)*(len_3p));
      if (var->flank3p==NULL) 
	{
	  free(var->flank5p);
	  free(var->one_allele);
	  free(var->other_allele);
	  free(var);
	  return NULL;
	}

      
      var->flank5p_or    = (Orientation*) malloc(sizeof(Orientation) * len_5p);
      if (var->flank5p_or==NULL)
	{
	  free(var->flank5p);
	  free(var->flank3p);
	  free(var->one_allele);
	  free(var->other_allele);
	  free(var);
	  return NULL;
	}
       var->one_allele_or = (Orientation*) malloc(sizeof(Orientation) * len_br1);
      if (var->one_allele_or==NULL)
	{
	  free(var->flank5p);
	  free(var->flank3p);
	  free(var->one_allele);
	  free(var->other_allele);
	  free(var->flank5p_or);
	  free(var);
	  return NULL;
	}
      var->other_allele_or = (Orientation*) malloc(sizeof(Orientation) * len_br2);
      if (var->other_allele_or==NULL)
	{
	  free(var->flank5p);
	  free(var->flank3p);
	  free(var->one_allele);
	  free(var->other_allele);
	  free(var->flank5p_or);
	  free(var->one_allele_or);
	  free(var);
	  return NULL;
	}
      var->flank3p_or    = (Orientation*) malloc(sizeof(Orientation) * len_3p);
      if (var->flank3p_or==NULL)
	{
	  free(var->flank5p);
	  free(var->flank3p);
	  free(var->one_allele);
	  free(var->other_allele);
	  free(var->flank5p_or);
	  free(var->one_allele_or);
	  free(var->other_allele_or);
	  free(var);
	  return NULL;
	}
      
      var->var_name=(char*) malloc(sizeof(char)*MAX_VARNAME_LEN);
      if (var->var_name==NULL)
	{
	  free(var->flank5p);
	  free(var->flank3p);
	  free(var->one_allele);
	  free(var->other_allele);	
	  free(var->flank5p_or);
	  free(var->one_allele_or);
	  free(var->other_allele_or);
	  free(var->flank3p_or);
	  free(var);
	  return NULL;
	}

      var->seq5p=(char*) malloc(sizeof(char)*(len_5p+kmer_size));
      if (var->seq5p==NULL)
	{
	  free(var->flank5p);
	  free(var->flank3p);
	  free(var->one_allele);
	  free(var->other_allele);	
	  free(var->flank5p_or);
	  free(var->one_allele_or);
	  free(var->other_allele_or);
	  free(var->flank3p_or);
	  free(var->var_name);
	  free(var);	  
	  return NULL;

	}

      var->seq_one=(char*) malloc(sizeof(char)*len_br1);
      if (var->seq_one==NULL)
	{
	  free(var->flank5p);
	  free(var->flank3p);
	  free(var->one_allele);
	  free(var->other_allele);	
	  free(var->flank5p_or);
	  free(var->one_allele_or);
	  free(var->other_allele_or);
	  free(var->flank3p_or);
	  free(var->var_name);
	  free(var->seq5p);
	  free(var);	  
	  return NULL;

	}


      var->seq_other=(char*) malloc(sizeof(char)*len_br2);
      if (var->seq_other==NULL)
	{
	  free(var->flank5p);
	  free(var->flank3p);
	  free(var->one_allele);
	  free(var->other_allele);	
	  free(var->flank5p_or);
	  free(var->one_allele_or);
	  free(var->other_allele_or);
	  free(var->flank3p_or);
	  free(var->var_name);
	  free(var->seq5p);
	  free(var->seq_one);
	  free(var);	  
	  return NULL;

	}

      var->seq3p =(char*) malloc(sizeof(char)*len_3p);
      if (var->seq3p==NULL)
	{
	  free(var->flank5p);
	  free(var->flank3p);
	  free(var->one_allele);
	  free(var->other_allele);	
	  free(var->flank5p_or);
	  free(var->one_allele_or);
	  free(var->other_allele_or);
	  free(var->flank3p_or);
	  free(var->var_name);
	  free(var->seq5p);
	  free(var->seq_one);
	  free(var->seq_other);
	  free(var);	  
	  return NULL;

	}


      var->len_flank5p=0;
      var->len_flank3p=0;
      var->len_one_allele=0;
      var->len_other_allele=0;
      var->which = unknown;
      var->var_name[0]='\0';
      var->seq5p[0]='\0';
      var->seq_one[0]='\0';
      var->seq_other[0]='\0';
      var->seq3p[0]='\0';

      return var;
    }
}

void free_VariantBranchesAndFlanks_object(VariantBranchesAndFlanks* var)
{
  free(var->flank5p);
  free(var->flank3p);
  free(var->one_allele);
  free(var->other_allele);
  free(var->flank5p_or);
  free(var->flank3p_or);
  free(var->one_allele_or);
  free(var->other_allele_or);
  free(var->seq5p);
  free(var->seq_one);
  free(var->seq_other);
  free(var->seq3p);
  free(var->var_name);
  free(var);
}


//use this for a var you have created on the stack
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

// suppose branch1 was node1 node2 node3. normally yoiu give an aray of node*'s  starting at node1
// this function allows you to pass in the array (node3, node2, node1) + the info that it is in reverse order
void set_alloced_variant_branches_and_flanks_allowing_inputargs_in_either_order(VariantBranchesAndFlanks* var, 
										dBNode** flank5p,    Orientation* flank5p_or,    int len_flank5p,  Orientation arraydir_5p, 
										dBNode** one_allele, Orientation* one_allele_or, int len_one_allele, Orientation arraydir_one,
										dBNode** other_allele, Orientation* other_allele_or, int len_other_allele, Orientation arraydir_other,
										dBNode** flank3p,    Orientation* flank3p_or,    int len_flank3p, Orientation arraydir_3p,
										WhichAlleleIsRef which)
{
  int i,j;
  if (arraydir_5p==reverse)
    {
      j=0;
      for (i=len_flank5p; i>=0; i--)
	{
	  var->flank5p[j]   =flank5p[i];
	  var->flank5p_or[j]=flank5p_or[i];
	  j++;
	}
    }
  else
    {
      for (i=0; i<=len_flank5p; i++)
	{
	  var->flank5p[i]   =flank5p[i];
          var->flank5p_or[i]=flank5p_or[i];
	}
    }
  var->len_flank5p   = len_flank5p;

  
  if (arraydir_one==reverse)
    {
      j=0;
      for (i=len_one_allele; i>=0; i--)
	{
	  var->one_allele[j]    = one_allele[i];
	  var->one_allele_or[j] = one_allele_or[i];
	  j++;
	}
    }
  else
    {
      for (i=0; i<=len_one_allele; i++)
	{
	  var->one_allele[i]    = one_allele[i];
	  var->one_allele_or[i] = one_allele_or[i];	  
	}
    }
  var->len_one_allele= len_one_allele;

  if (arraydir_other==reverse)
    {
      j=0;
      for (i=len_other_allele; i>=0; i--)
	{
	  var->other_allele[j]    = other_allele[i];
	  var->other_allele_or[j] = other_allele_or[i];
	  j++;
	}
    }
  else
    {
      for (i=len_other_allele; i>=0; i--)
	{
	  var->other_allele[i]    = other_allele[i];
	  var->other_allele_or[i] = other_allele_or[i];
	}      
    }
  var->len_other_allele= len_other_allele;

  if (arraydir_3p==reverse)
    {
      j=0;
      for (i=len_flank3p; i>=0; i--)
	{
	  var->flank3p[j]       =flank3p[i];
	  var->flank3p_or[j]    =flank3p_or[i];
	  j++;
	}
    }
  else
    {
      for (i=len_flank3p; i>=0; i--)
        {
          var->flank3p[i]       =flank3p[i];
          var->flank3p_or[i]    =flank3p_or[i];
	}
    }
  var->len_flank3p   = len_flank3p;
  var->which         = which;
}

void set_status_of_nodes_in_branches(VariantBranchesAndFlanks* var, NodeStatus status)
{
  int i;

  for (i=0; i<var->len_one_allele; i++)
    {
      db_node_set_status( (var->one_allele)[i], status);
    }

  for (i=0; i<var->len_other_allele; i++)
    {
      db_node_set_status( (var->other_allele)[i], status);
    }

}


void set_status_of_genotyping_nodes_in_branches(GenotypingVariantBranchesAndFlanks* var, NodeStatus status)
{
  int i;

  for (i=0; i<var->len_one_allele; i++)
    {
      db_genotyping_node_set_status( (var->one_allele)[i], status);
    }

  for (i=0; i<var->len_other_allele; i++)
    {
      db_genotyping_node_set_status( (var->other_allele)[i], status);
    }

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

void set_variant_branches_but_flanks_to_null(VariantBranchesAndFlanks* var, 
					     dBNode** one_allele, Orientation* one_allele_or, int len_one_allele, 
					     dBNode** other_allele, Orientation* other_allele_or, int len_other_allele, 
					     WhichAlleleIsRef which)
{
  var->flank5p       = NULL;
  var->flank5p_or    = NULL;
  var->len_flank5p   = 0;
  var->one_allele    = one_allele;
  var->one_allele_or = one_allele_or;
  var->len_one_allele= len_one_allele;
  var->other_allele    = other_allele;
  var->other_allele_or = other_allele_or;
  var->len_other_allele= len_other_allele;
  var->flank3p       = NULL;
  var->flank3p_or    = NULL;
  var->len_flank3p   = 0;
  var->which         = which;
}


void set_genotyping_variant_branches_but_flanks_to_null(GenotypingVariantBranchesAndFlanks* var, 
							GenotypingElement** one_allele, Orientation* one_allele_or, int len_one_allele, 
							GenotypingElement** other_allele, Orientation* other_allele_or, int len_other_allele, 
							WhichAlleleIsRef which)
{
  var->flank5p       = NULL;
  var->flank5p_or    = NULL;
  var->len_flank5p   = 0;
  var->one_allele    = one_allele;
  var->one_allele_or = one_allele_or;
  var->len_one_allele= len_one_allele;
  var->other_allele    = other_allele;
  var->other_allele_or = other_allele_or;
  var->len_other_allele= len_other_allele;
  var->flank3p       = NULL;
  var->flank3p_or    = NULL;
  var->len_flank3p   = 0;
  var->which         = which;
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


void get_all_genotype_log_likelihoods_at_bubble_call_for_one_colour(AnnotatedPutativeVariant* annovar, double seq_error_rate_per_base, 
								    double sequencing_depth_of_coverage, int read_length, int colour)
{
  boolean too_short = false;
  int initial_covg_plus_upward_jumps_branch1 = 
    count_reads_on_allele_in_specific_colour(annovar->var->one_allele, annovar->var->len_one_allele, colour, &too_short);
  int initial_covg_plus_upward_jumps_branch2 = 
    count_reads_on_allele_in_specific_colour(annovar->var->other_allele, annovar->var->len_other_allele, colour, &too_short);

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

  double theta_one   = ((double)(sequencing_depth_of_coverage * annovar->var->len_one_allele))  /( (double) read_length );
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






void get_all_haploid_genotype_log_likelihoods_at_bubble_call_for_one_colour(AnnotatedPutativeVariant* annovar, double seq_error_rate_per_base, 
									    double sequencing_depth_of_coverage, int read_length, int colour)
{
  boolean too_short = false;
  int initial_covg_plus_upward_jumps_branch1 = 
    count_reads_on_allele_in_specific_colour(annovar->var->one_allele, annovar->var->len_one_allele, colour, &too_short);
  int initial_covg_plus_upward_jumps_branch2 = 
    count_reads_on_allele_in_specific_colour(annovar->var->other_allele, annovar->var->len_other_allele, colour, &too_short);

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
  annovar->gen_log_lh[colour].log_lh[het]       = annovar->gen_log_lh[colour].log_lh[hom_one] + annovar->gen_log_lh[colour].log_lh[hom_other] - 99999;



}

/*
void get_all_haploid_genotype_log_likelihoods_at_non_SNP_PD_call_for_one_colour(AnnotatedPutativeVariant* annovar, double seq_error_rate_per_base, double sequencing_depth_of_coverage, int read_length, int colour)
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
  //if not a SNP
  else if ( (annovar->var->len_one_allele != annovar->var->len_other_allele) || (annovar->var->len_one_allele>annovar->kmer+1) )
    {
      //set to hom_other (the variant allele)
      annovar->gen_log_lh[colour].log_lh[hom_one]=-999;
      annovar->gen_log_lh[colour].log_lh[het]=-9999999;
      annovar->gen_log_lh[colour].log_lh[hom_other]=0;
      return;
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
  annovar->gen_log_lh[colour].log_lh[het] = annovar->gen_log_lh[colour].log_lh[hom_one] + annovar->gen_log_lh[colour].log_lh[hom_other] - 999999;



}
*/

//assuming a pair of branches really do make up a variant, calculate the log likelihood of a genotype
//under the model described in our paper (eg used for HLA)
//theta here is as used in the paper: (D/R) * length of branch/allele. NOT the same theta as seen in model_selection.c
//assumes called by BubbleCaller, so no overlaps between alleles.
double get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(zygosity genotype, double error_rate_per_base, int covg_branch_1, int covg_branch_2, 
									double theta_one, double theta_other, int kmer)
{
  //  printf("DEBUG covg branch 1 is %d, and branch2 is %d, and theta 1,2 are %f, %f, and error ratew per base is %f\n", covg_branch_1, covg_branch_2, theta_one, theta_other, error_rate_per_base);
  if (genotype==hom_one)
    {
      //Apply formula for likelihood in section 9.0 of Supp. Methods of paper; no unique segment, one shared segment
      return covg_branch_1*log(theta_one) - theta_one
	- gsl_sf_lnfact(covg_branch_1)  + covg_branch_2 *log(error_rate_per_base) ;
    }
  else if (genotype==hom_other)
    {
      //Apply formula for likelihood in section 9.0 of Supp. Methods of paper; no unique segment, one shared segment
      return covg_branch_2*log(theta_other) - theta_other
	- gsl_sf_lnfact(covg_branch_2) + covg_branch_1 *log(error_rate_per_base)  ;
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
  gl->log_lh[hom_one]=0;
  gl->log_lh[het]=0;
  gl->log_lh[hom_other]=0;
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
//note I do not want to create an array on the stack - these things can be very long
// so relies on the prealloced array of dBNode* 's passed in
// annoying that can't use templates or something - see below for a similar function with different input
int count_reads_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short)
{

  if ( (len==0) || (len==1) )
    {
      *too_short=true;
      return 0;//previously this returned -1, which seems perverse to me
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



//WARNING - this is for use when we dissect an allele into subchunks, so here we do not want to be ignoring first and last elements (cf above)
int count_reads_on_allele_in_specific_colour_given_array_of_cvgs(int* covgs, int len, boolean* too_short)
{

  if ( (len==0) || (len==1) )
    {
      *too_short=true;
      return -1;
    }


  int total= covgs[0];

  int i;

  for (i=1; i<len; i++)
    {
      int jump = covgs[i] - covgs[i-1];

      //we add a little check to ensure that we ignore isolated nodes with higher covg - we count jumps only if they are signifiers of a new read arriving
      // and one node does not a read make
      int diff_between_next_and_prev=-1;
      if (i<len-1)
	{
	  diff_between_next_and_prev=  covgs[i+1]-covgs[i-1];
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

