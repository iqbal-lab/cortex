#include <db_variants.h>


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
