#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <dB_graph.h>

//for fasta


extern int MAX_FILENAME_LENGTH;
extern int MAX_READ_LENGTH;

int load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * count_kmers, long long * bad_reads, int max_read_length, dBGraph* db_graph, EdgeArrayType type, 
								       int index);
int load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long * count_kmers, long long * bad_reads,  char quality_cut_off, int max_read_length, dBGraph* db_graph,
								       EdgeArrayType type, int index);

long long load_population_as_fasta(char* filename,  long long* count_kmers, long long* bad_reads, dBGraph* db_graph);
int load_all_fasta_for_given_person_given_filename_of_file_listing_their_fasta_files(char* f_name, long long* count_kmers, long long* bad_reads, dBGraph* db_graph, int index);

long long load_population_as_fastq(char* filename,  long long* count_kmers, long long* bad_reads, char quality_cutoff, dBGraph* db_graph);
int load_all_fastq_for_given_person_given_filename_of_file_listing_their_fastq_files(char* f_name, long long* count_kmers, long long* bad_reads, char quality_cutoff, dBGraph* db_graph, int index);






#endif /* FILE_READER_H_ */
