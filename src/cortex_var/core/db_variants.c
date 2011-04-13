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

#include <db_variants.h>
#include <dB_graph.h>
#include <dB_graph_population.h>


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
