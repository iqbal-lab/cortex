#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <dB_graph.h>

int load_fasta_data_from_filename_into_graph(char* filename, long long *count_kmers, char quality_cut_off, int max_read_length, dBGraph* db_graph);

//for fastq
int load_fastq_data_from_filename_into_graph(char* filename, long long * count_kmers, char quality_cut_off, int max_read_length, dBGraph* db_graph);

#endif /* FILE_READER_H_ */
