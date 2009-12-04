#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <dB_graph.h>
#include <seq.h>

int MAX_READ_LENGTH;

extern int MAX_FILENAME_LENGTH;


//these routines return the length of the read sequence, for the binary file is all the kmers conctenated

//for short fasta entries - reads or similar 

long long load_fasta_from_filename_into_graph(char* filename, long long * bad_reads, long long* dup_reads, int max_chunk_length, 
					      boolean remove_duplicates_single_endedly, dBGraph* db_graph);

//for fastq
long long load_fastq_from_filename_into_graph(char* filename, long long * bad_reads,  char quality_cut_off, long long* dup_reads, int max_read_length, 
					      boolean remove_duplicates_single_endedly, dBGraph* db_graph);

long long load_paired_fastq_from_filenames_into_graph(char* filename1, char* filename2, long long * bad_reads,  char quality_cut_off, int max_read_length,
                                                      long long* dup_reads, boolean remove_duplicates, dBGraph* db_graph);
long long load_list_of_paired_end_fastq_into_graph(char* list_of_left_mates, char* list_of_right_mates, char quality_cut_off, int max_read_length,
                                                   long long* bad_reads, long long* num_dups, int* count_file_pairs, boolean remove_dups, dBGraph* db_graph);

//for reference

void read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
                                                                            long long * bad_reads, char quality_cut_off, int max_read_length, dBGraph * db_graph);

void read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(char* f_name, dBGraph* db_graph);

void read_all_ref_chromosomes_and_mark_graph(dBGraph* db_graph);



// This would be in seq.c, but that has no knowledge of dBGraph
int get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(char * seq,  char * qualities, int length, char quality_cut_off,
                                                                                  KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph);



int get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(char * seq,  char * qualities, int length, char quality_cut_off,
                                                                                     KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph);




// for dumping clean fasta files from fastq - ie only holding reads that lie entirely in a (presumably cleaned) graph
void read_fastq_and_print_reads_that_lie_in_graph(FILE* fp, FILE* fout, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
						  long long * bad_reads, int max_read_length, dBGraph * db_graph,
						  boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);

void read_fastq_and_print_subreads_that_lie_in_graph_breaking_at_edges_or_kmers_not_in_graph(FILE* fp, FILE* fout,
                                                                                             int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry,
                                                                                                                 boolean * full_entry),
                                                                                             long long * bad_reads, int max_read_length, dBGraph * db_graph,
                                                                                             boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);


//for binary
long long  load_binary_from_filename_into_graph(char* filename, dBGraph* db_graph);




#endif /* FILE_READER_H_ */
