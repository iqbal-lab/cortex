#include <element.h>


#ifndef DB_VARIANTS_H_
#define DB_VARIANTS_H_

//a variant can be hom for oen allele or the other, het or absent in an individual
typedef enum
  {
    hom_one = 0,
    hom_other = 1,
    het = 2,
    absent =3
  }zygosity;

typedef enum
 {
    full_flank_format   = 0,
    glf  = 1,
 } Variant_File_Format ;

typedef enum
  {
    allele_one = 0,
    allele_other = 1,
  } WhichAllele;

typedef enum
  {
    first = 0,
    unknown = 1,
  } WhichAlleleIsRef ;

typedef struct{
  dBNode** flank5p;
  Orientation* flank5p_or;
  int len_flank5p;
  dBNode** one_allele; 
  Orientation* one_allele_or; 
  int len_one_allele;
  dBNode** other_allele; 
  Orientation* other_allele_or;
  int len_other_allele;
  dBNode** flank3p;    
  Orientation* flank3p_or;
  int len_flank3p;
  WhichAlleleIsRef which;
} VariantBranchesAndFlanks;


void set_variant_branches_and_flanks(VariantBranchesAndFlanks* var, 
				     dBNode** flank5p,    Orientation* flank5p_or,    int len_flank5p,
				     dBNode** one_allele, Orientation* one_allele_or, int len_one_allele, 
				     dBNode** other_allele, Orientation* other_allele_or, int len_other_allele, 
				     dBNode** flank3p,    Orientation* flank3p_or,    int len_flank3p, WhichAlleleIsRef which);


void exact_copy_variant_branches_and_flanks(VariantBranchesAndFlanks copy_to, const VariantBranchesAndFlanks copy_from);
void copy_variant_branches_and_flanks_switching_branches(VariantBranchesAndFlanks copy_to, const VariantBranchesAndFlanks copy_from);
void action_set_flanks_and_branches_to_be_ignored(VariantBranchesAndFlanks* var);
void db_variant_action_do_nothing(VariantBranchesAndFlanks* var);

boolean  db_variant_precisely_one_allele_is_in_given_func_of_colours(VariantBranchesAndFlanks* var, Edges (*get_colour)(const dBNode*), dBGraph* db_graph, WhichAllele* which,
                                                                     int num_bases_agreement_in_flank_required );

zygosity db_variant_get_zygosity_in_given_func_of_colours(VariantBranchesAndFlanks* var, Edges (*get_colour)(const dBNode*), dBGraph* db_graph);
#endif
