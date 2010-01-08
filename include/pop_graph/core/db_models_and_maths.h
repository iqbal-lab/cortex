#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <seq.h>
#include <string.h>
#include <limits.h>
#include <file_reader.h>
#include <stdlib.h>
#include <math.h>


#ifndef DB_MODELS_AND_MATHS_
#define DB_MODELS_AND_MATHS_


extern const double LARGE_NEGATIVE_NUMBER;

int factorial( int n );
int max(int a, int b, int c, int d);


// Given a sequence of coverage values for a set of contiguous nodes,
// where we have factored out effects due to rest of genome or repeat/multiplicity
//calculate probability, under a Poisson approx to a binomial model, 
// that we see those coverages
double get_probability_of_observed_covgs(int* coverages, int array_len, double lambda);
//same for log prob - expect to have less trunc/rounding errors
double get_log_probability_of_observed_covgs(int* coverages, int array_len, double lambda);

void get_normalised_coverage(int population_size, int avg_depth_of_covg, int avg_read_len, 
			    int* ref_normalised_covg, int* alt_normalised_covg,
			    int* ref_multiplicity_in_ref, int* alt_multiplicity_in_ref,
			    VariantBranchesAndFlanks* var, dBGraph* db_graph, int colour, int colour_ref);

double calc_log_likelihood_of_data_seen_on_one_allele_excluding_nodes_on_both_alleles(int* array_of_covgs, int* array_of_mult_wrt_self, int* array_of_mult_wrt_other, int* normalised_covg,
										      int len_arrays, int allele_freq, double avg_depth_of_covg_per_haploid, int avg_read_len, short kmer_size);

//use this when you're considering ref-ref or alt-alt models, so there is no "other" allele
double calc_log_likelihood_of_data_seen_on_one_allele(int* array_of_covgs, int* array_of_mult_wrt_self,  int* normalised_covg, 
						      int len_arrays, int allele_freq, double avg_depth_of_covg_per_haploid, int avg_read_len, short kmer_size);




#endif
