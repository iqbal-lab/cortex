#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <seq.h>
#include <dB_graph.h>
#include <pop_globals.h>

//for fasta


extern int MAX_FILENAME_LENGTH;
extern int MAX_READ_LENGTH;

int load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * bad_reads, int max_read_length, dBGraph* db_graph, EdgeArrayType type, 
								       int index);
int load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * bad_reads,  char quality_cut_off, int max_read_length, dBGraph* db_graph,
								       EdgeArrayType type, int index);


long long load_population_as_fasta(char* filename, long long* bad_reads, dBGraph* db_graph);
int load_all_fasta_for_given_person_given_filename_of_file_listing_their_fasta_files(char* f_name, long long* bad_reads, dBGraph* db_graph, int index);

long long load_population_as_fastq(char* filename, long long* bad_reads, char quality_cutoff, dBGraph* db_graph);
int load_all_fastq_for_given_person_given_filename_of_file_listing_their_fastq_files(char* f_name, long long* bad_reads, char quality_cutoff, dBGraph* db_graph, int index);

int load_chromosome_overlap_data(char* f_name,  dBGraph* db_graph, int which_chromosome);

int load_sv_trio_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph);
int load_individual_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph, EdgeArrayType type, int index);

//gets number_of_bases_to_load's worth of kmers, and returns the corresponding nodes, orientations etc in he array passed in.

int load_seq_into_array(FILE* chrom_fptr, int number_of_bases_to_load, int length_of_arrays,
			    dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
			    Sequence* seq, KmerSlidingWindow* kmer_window, boolean expecting_new_fasta_entry,  dBGraph * db_graph);


#endif /* FILE_READER_H_ */
