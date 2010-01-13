#include <db_models_and_maths.h>

#ifndef DB_GENOTYPING_H_
#define DB_GENOTYPING_H_





void get_node_multiplicities(dBNode** branch1, int len_branch1, dBNode** branch2, int len_branch2,
                             int** br1_multiplicity_in_br1, int** br2_multiplicity_in_br2,
                             int** br1_multiplicity_in_br2, int** br2_multiplicity_in_br1,
			     boolean only_count_nodes_with_edge_in_specified_colour_func,
                             Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*) );

void genotype_all_variants_in_fff_file(char* filename, Variant_File_Format output_format, char* out_filename, int max_branch_len,
				       void (*genotype_model) (VariantBranchesAndFlanks*, dBGraph*, int, Variant_File_Format, FILE*, 
							       int, int, int, double), 
				       dBGraph* db_graph, int colour1 , int colour_ref, int population_size_colour1, 
				       double avg_depth_of_covg_pop1, int avg_read_len);

void genotype_model_simple(VariantBranchesAndFlanks* var, dBGraph* db_graph, int avg_read_len,
			   Variant_File_Format output_format, FILE* output_fptr, 
			   int colour1 , int colour_ref, int population_size_colour1, double avg_depth_of_covg_pop1);


#endif
