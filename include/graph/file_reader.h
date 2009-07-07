#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <dB_graph.h>
#include <seq.h>

int MAX_READ_LENGTH;


//these routines return the length of the read sequence, for the binary file is all the kmers conctenated

//for short fasta entries - reads or similar 

long long load_fasta_from_filename_into_graph(char* filename, long long * bad_reads, int max_chunk_length, dBGraph* db_graph);

//for fastq
long long load_fastq_from_filename_into_graph(char* filename, long long * bad_reads,  char quality_cut_off, int max_read_length, dBGraph* db_graph);

//for reference

void read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
                                                                            long long * bad_reads, char quality_cut_off, int max_read_length, dBGraph * db_graph);

void read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(char* f_name, dBGraph* db_graph);

void read_all_ref_chromosomes_and_mark_graph(dBGraph* db_graph);




//for binary
long long  load_binary_from_filename_into_graph(char* filename, dBGraph* db_graph);




#endif /* FILE_READER_H_ */
