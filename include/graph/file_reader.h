#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <dB_graph.h>
#include <seq.h>

//these routines return the length of the read sequence, for the binary file is all the kmers conctenated

//for fasta
long long load_fasta_data_from_filename_into_graph(char* filename, long long * bad_reads, int max_read_length, dBGraph* db_graph);

//for fastq
long long load_fastq_data_from_filename_into_graph(char* filename, long long * bad_reads,  char quality_cut_off, int max_read_length, dBGraph* db_graph);

//for binary
long long  load_binary_data_from_filename_into_graph(char* filename, dBGraph* db_graph);

long long load_binary_data_from_filename_into_graph_with_unique_kmers(char* filename, dBGraph* db_graph);

//for reference
void read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length),
                                                                            long long * bad_reads, char quality_cut_off, int max_read_length, dBGraph * db_graph);

void read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(char* f_name, dBGraph* db_graph);

void read_all_ref_chromosomes_and_mark_graph(dBGraph* db_graph);

#endif /* FILE_READER_H_ */
